/**********************************************************************************************
© 2020. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and
to permit others to do so.
This program is open source under the BSD-3 License.
**********************************************************************************************/

/// @file   logger.hpp
/// @brief  MPI- and GPU-safe printf-style logger for Fierro.
///
/// ============================================================================
/// WHAT THIS FILE GIVES YOU
/// ============================================================================
///
///   fierro::Logger        -- owning, non-copyable, one-per-Driver log sink.
///                            Host-only. Not capturable into kernels.
///   fierro::Logger::Handle -- trivially-copyable POD view of a Logger, safe to
///                            capture by value into FOR_ALL / KOKKOS_LAMBDA on
///                            every Kokkos backend (Serial/OpenMP/CUDA/HIP/SYCL).
///   fierro::log_handle(Logger*) -- null-safe shorthand for `p ? p->handle() :
///                            Handle{}`; use with Solver::log_handle() as
///                            `const auto log_dev = log_handle();` before FOR_ALL.
///   FLOG_DEV(handle, LEVEL, fmt, ...) -- optional single-printf macro for
///                            performance-critical device log sites.
///
/// ============================================================================
/// RULES (read before using)
/// ============================================================================
///
///  1) Host side: call   log.info("rank-local msg %d\n", x);
///                anywhere on the host, including inside loops whose iteration
///                count varies across MPI ranks. These calls do NOT touch MPI.
///
///  2) Device side: obtain a POD view with   auto lh = log.handle();
///                then capture `lh` by value into your FOR_ALL / KOKKOS_LAMBDA
///                and call lh.info(...), lh.warn(...), etc. Never try to
///                capture the owning Logger itself -- it is non-copyable.
///
///  3) log.flush() is the ONLY collective operation in this class. Call it at
///                a program point every rank reaches (end of a phase, end of a
///                time step). NEVER call flush inside a loop whose trip count
///                varies across ranks.
///
///  4) Never call any MPI_* from inside a FOR_ALL / KOKKOS_LAMBDA.
///
///  5) log.set_level(...) / log.set_file_prefix(...) take effect for future
///     logging only. They do NOT retroactively change Handles already captured
///     into running kernels -- the Handle carries a snapshot of `level`.
///
///  6) flush() is NOT host-thread-safe. The hot-path append IS (it uses
///     Kokkos::atomic_fetch_add for serialization). Call flush from a single
///     host thread.
///
/// ============================================================================
/// PORTABILITY
/// ============================================================================
///
///   The Logger's state uses only Kokkos-portable primitives:
///     - Kokkos::View<char*,   Kokkos::HostSpace>  (the rank-local buffer)
///     - Kokkos::View<size_t,  Kokkos::HostSpace>  (atomic byte counter)
///     - Kokkos::View<size_t,  Kokkos::HostSpace>  (atomic drop counter)
///     - Kokkos::atomic_fetch_add                  (replaces std::mutex)
///     - raw char*, size_t*, int, LogLevel, FILE*, MPI_Comm        (POD)
///
///   No std::string, std::vector, std::mutex, std::thread, or std::shared_ptr
///   appears anywhere in the class state or in flush(). All flush()
///   temporaries are Kokkos::View<..., Kokkos::HostSpace>.
///
/// ============================================================================
/// OUTPUT FORMAT (both stdout and per-rank files)
/// ============================================================================
///
///   Every line is tagged with rank and 5-char padded level:
///
///       [rank 0][INFO ] Starting SGH3D setup
///       [rank 0][TRACE] entering calc_corner_mass
///       [rank 2][WARN ] 3 elements near zero volume in material 2
///       [rank 1][ERROR] negative density detected in elem 4211
///
/// ============================================================================
/// FILE OUTPUT
/// ============================================================================
///
///   Each rank opens its own file at construction: "<file_prefix><rank>" with
///   mode "w" (truncating). Default file_prefix is "Fierro_log_", yielding
///   Fierro_log_0, Fierro_log_1, ... . Change via the ctor argument or
///   set_file_prefix(). Pass nullptr to disable file output.
///
///   File I/O in flush() uses std::fwrite on the rank's own FILE*, so it never
///   involves MPI and is unaffected by rank-varying loop bounds.

#ifndef FIERRO_LOGGER_H
#define FIERRO_LOGGER_H

#include <Kokkos_Core.hpp>
#include <Kokkos_Printf.hpp>
#include <mpi.h>

#include <cstdio>    // std::snprintf, std::fopen, std::fwrite, FILE*
#include <cstring>   // std::memcpy
#include <cstddef>   // size_t
#include <type_traits>

/// @def   FIERRO_LOG_DEVICE_MAX_LEVEL
/// @brief Minimum level that actually emits a Kokkos::printf on the device.
///        Levels below this are stripped via `if constexpr` in the level-sugar
///        methods (info/warn/error/debug/trace), producing zero device-side
///        instructions at call sites that are compiled out.
///
/// Defaults: release (NDEBUG) -> Warn, debug -> Trace.
/// Override on the command line:
///     -DFIERRO_LOG_DEVICE_MAX_LEVEL=::fierro::LogLevel::Info
#ifndef FIERRO_LOG_DEVICE_MAX_LEVEL
#  ifdef NDEBUG
#    define FIERRO_LOG_DEVICE_MAX_LEVEL ::fierro::LogLevel::Warn
#  else
#    define FIERRO_LOG_DEVICE_MAX_LEVEL ::fierro::LogLevel::Trace
#  endif
#endif

namespace fierro {

/// @brief Severity level ordering: Trace < Debug < Info < Warn < Error < Off.
///
/// Host-side filtering: a call with `lvl < Logger::level()` returns early
/// before any snprintf. Device-side filtering: the level-sugar methods
/// (info/warn/error/debug/trace) compile out below FIERRO_LOG_DEVICE_MAX_LEVEL
/// via `if constexpr`, so disabled levels emit zero instructions in
/// CUDA/HIP/SYCL kernels. Off disables all logging on the host.
enum class LogLevel : int {
    Trace = 0,
    Debug = 1,
    Info  = 2,
    Warn  = 3,
    Error = 4,
    Off   = 5
};

/// @brief Pad a LogLevel to a 5-char tag so "[INFO ]", "[WARN ]", "[ERROR]",
///        "[TRACE]", "[DEBUG]" visually align. KOKKOS_INLINE_FUNCTION so it is
///        callable on host and device at zero cost (resolves at compile time
///        for literal arguments).
KOKKOS_INLINE_FUNCTION
constexpr const char* level_tag(LogLevel l) {
    switch (l) {
        case LogLevel::Trace: return "TRACE";
        case LogLevel::Debug: return "DEBUG";
        case LogLevel::Info : return "INFO ";
        case LogLevel::Warn : return "WARN ";
        case LogLevel::Error: return "ERROR";
        case LogLevel::Off  : return "OFF  ";
    }
    return "?????";
}


// ============================================================================
//  Logger  (declaration + inline definitions)
// ============================================================================

class Logger {
public:
    // ------------------------------------------------------------------------
    //  Handle -- trivially-copyable POD view suitable for kernel capture.
    // ------------------------------------------------------------------------

    /// @brief Lightweight, trivially-copyable view of a Logger for capture by
    ///        value into FOR_ALL / KOKKOS_LAMBDA on every Kokkos backend.
    ///
    /// All members are POD pointers or ints. Never touches std::string/vector/
    /// mutex. Obtain via Logger::handle() on the host; then capture a local
    /// copy into the kernel. Never try to capture the owning Logger itself --
    /// it is non-copyable.
    ///
    /// @warning The `level` field is a snapshot of Logger::level() at the
    ///          moment Handle was produced. Changing Logger::set_level(...)
    ///          after a Handle has been captured into a running kernel does
    ///          not affect that kernel's filtering. Refresh with
    ///          `log.handle()` after any set_level call.
    struct Handle {
        char*    buf_ptr  = nullptr;              ///< Host buffer, never deref'd on device
        size_t*  used     = nullptr;              ///< Atomic byte counter (host)
        size_t*  dropped  = nullptr;              ///< Atomic overflow counter (host)
        size_t   cap      = 0;                    ///< Buffer capacity in bytes
        int      rank     = 0;                    ///< MPI rank in Logger::comm()
        LogLevel level    = LogLevel::Info;       ///< Minimum level (snapshot)

        /// @brief Append one printf-style log line. Unified host/device.
        ///
        /// Host path: formats "[rank N][LEVEL] <user message>" into a 512-byte
        /// stack buffer, then reserves space in the rank-local ring buffer
        /// with Kokkos::atomic_fetch_add and memcpy's in. If the reservation
        /// exceeds `cap`, the message is dropped and the atomic drop counter
        /// is incremented (no rollback -- rollback races with other reservers).
        ///
        /// Device path: emits two Kokkos::printf calls (tag then body). For
        /// performance-critical device sites, prefer the fused FLOG_DEV macro.
        ///
        /// @tparam Args   Types of the printf arguments (deduced).
        /// @param  lvl    Severity of this message. Host filters against
        ///                this Handle's snapshot of `level`; device path is
        ///                unfiltered at runtime (use compile-time gating via
        ///                the sugar methods + FIERRO_LOG_DEVICE_MAX_LEVEL).
        /// @param  fmt    printf-compatible format string literal. Must
        ///                outlive this call.
        /// @param  args   printf-compatible arguments matching `fmt`.
        ///
        /// @note Non-collective. Does NOT call MPI.
        /// @note Host-thread-safe via Kokkos::atomic_fetch_add.
        /// @note If the message would push the buffer past `cap`, it is
        ///       dropped and a warning is emitted during the next flush().
        ///
        /// @code
        ///     // Host:
        ///     log.log(fierro::LogLevel::Info, "num_elems=%zu\n", n);
        ///
        ///     // Device:
        ///     auto lh = log.handle();
        ///     FOR_ALL(i, 0, N, {
        ///         lh.log(fierro::LogLevel::Warn, "i=%d\n", (int)i);
        ///     });
        /// @endcode
        template<typename... Args>
        KOKKOS_INLINE_FUNCTION
        void log(LogLevel lvl, const char* fmt, Args... args) const {
            KOKKOS_IF_ON_DEVICE((
                if (cap == 0u) return;
                // Generic device path: two Kokkos::printf calls.
                // Prefer FLOG_DEV for hot sites to avoid cross-thread
                // interleaving between the tag and the body.
                Kokkos::printf("[rank %d][%s] ", rank, level_tag(lvl));
                Kokkos::printf(fmt, args...);
            ))
            KOKKOS_IF_ON_HOST((
                if (static_cast<int>(lvl) < static_cast<int>(level)) return;
                append_host_(lvl, fmt, args...);
            ))
        }

        /// @brief Sugar for log(Trace, ...). On device, compiles out below
        ///        FIERRO_LOG_DEVICE_MAX_LEVEL via `if constexpr` (zero
        ///        instructions when disabled).
        template<typename... Args>
        KOKKOS_INLINE_FUNCTION
        void trace(const char* fmt, Args... args) const {
            leveled_<LogLevel::Trace>(fmt, args...);
        }

        /// @brief Sugar for log(Debug, ...). Compile-time-gated on device.
        template<typename... Args>
        KOKKOS_INLINE_FUNCTION
        void debug(const char* fmt, Args... args) const {
            leveled_<LogLevel::Debug>(fmt, args...);
        }

        /// @brief Sugar for log(Info, ...). Compile-time-gated on device.
        template<typename... Args>
        KOKKOS_INLINE_FUNCTION
        void info(const char* fmt, Args... args) const {
            leveled_<LogLevel::Info>(fmt, args...);
        }

        /// @brief Sugar for log(Warn, ...). Compile-time-gated on device.
        template<typename... Args>
        KOKKOS_INLINE_FUNCTION
        void warn(const char* fmt, Args... args) const {
            leveled_<LogLevel::Warn>(fmt, args...);
        }

        /// @brief Sugar for log(Error, ...). Compile-time-gated on device.
        template<typename... Args>
        KOKKOS_INLINE_FUNCTION
        void error(const char* fmt, Args... args) const {
            leveled_<LogLevel::Error>(fmt, args...);
        }

        /// @brief False for a default-constructed handle (e.g. no logger attached).
        ///        Use `if (log_dev)` before calling trace/info/… on a captured Handle.
        KOKKOS_INLINE_FUNCTION explicit operator bool() const noexcept {
            return buf_ptr != nullptr && used != nullptr && dropped != nullptr &&
                   cap != 0u;
        }

    private:
        /// @brief Compile-time-typed shim used by info/warn/error/debug/trace.
        ///        Enables zero-cost device stripping via `if constexpr` against
        ///        FIERRO_LOG_DEVICE_MAX_LEVEL. On host, still runtime-filtered
        ///        against Logger::level() via the generic log() path.
        template<LogLevel Lvl, typename... Args>
        KOKKOS_INLINE_FUNCTION
        void leveled_(const char* fmt, Args... args) const {
            KOKKOS_IF_ON_DEVICE((
                if (cap == 0u) return;
                if constexpr (static_cast<int>(Lvl) >=
                              static_cast<int>(FIERRO_LOG_DEVICE_MAX_LEVEL)) {
                    Kokkos::printf("[rank %d][%s] ", rank, level_tag(Lvl));
                    Kokkos::printf(fmt, args...);
                }
            ))
            KOKKOS_IF_ON_HOST((
                if (static_cast<int>(Lvl) < static_cast<int>(level)) return;
                append_host_(Lvl, fmt, args...);
            ))
        }

        /// @brief Host-only append: format once, reserve atomically, memcpy.
        ///        Overflow is counted, not stored (no rollback).
        template<typename... Args>
        inline void append_host_(LogLevel lvl, const char* fmt, Args... args) const {
            KOKKOS_IF_ON_HOST((
                if (buf_ptr == nullptr || used == nullptr || dropped == nullptr ||
                    cap == 0u) {
                    return;
                }
                char stack_buf[512];
                int pfx = std::snprintf(stack_buf, sizeof(stack_buf),
                                        "[rank %d][%s] ", rank, level_tag(lvl));
                if (pfx <= 0) return;
                int msg = std::snprintf(stack_buf + pfx,
                                        sizeof(stack_buf) - static_cast<size_t>(pfx),
                                        fmt, args...);
                if (msg <= 0) return;
                size_t want = static_cast<size_t>(pfx) + static_cast<size_t>(msg);
                if (want > sizeof(stack_buf)) want = sizeof(stack_buf);

                // Monotonically grow `used` via atomic reservation. Writes
                // only happen when the reservation fits in `cap`; otherwise
                // the message is dropped (no rollback -- rollback races with
                // other reservers).
                size_t off = Kokkos::atomic_fetch_add(used, want);
                if (off + want <= cap) {
                    std::memcpy(buf_ptr + off, stack_buf, want);
                } else {
                    Kokkos::atomic_fetch_add(dropped, static_cast<size_t>(1));
                }
            ))
        }
    };  // struct Handle

    // Compile-time guarantee the Handle is a POD suitable for kernel capture.
    static_assert(std::is_trivially_copyable<Handle>::value,
                  "fierro::Logger::Handle must be trivially copyable to be "
                  "captured by value into FOR_ALL / KOKKOS_LAMBDA.");

    // ------------------------------------------------------------------------
    //  Logger construction / destruction
    // ------------------------------------------------------------------------

    /// @brief Construct a Logger bound to an MPI communicator and open a
    ///        per-rank log file.
    ///
    /// @param comm            MPI communicator used by flush(). NOT
    ///                        duplicated; the caller must ensure `comm`
    ///                        outlives the Logger. Defaults to
    ///                        MPI_COMM_WORLD. MPI must already be
    ///                        initialized (MPI_Init/MPI_Init_thread).
    /// @param level           Initial host-side minimum level. Messages below
    ///                        `level` are dropped before any formatting.
    ///                        Default LogLevel::Info.
    /// @param file_prefix     Filename prefix for this rank's on-disk log
    ///                        file. The file opened is "<file_prefix><rank>"
    ///                        (mode "w", truncating). Pass nullptr to skip
    ///                        file output entirely. Default "Fierro_log_".
    /// @param capacity_bytes  Size of the rank-local buffer, in bytes.
    ///                        Default 1 MiB. Messages past this capacity
    ///                        between flushes are dropped (and counted).
    ///
    /// @note Host-only. Every rank must construct with the same `comm`; no
    ///       MPI calls are issued inside the ctor itself.
    /// @note Allocates three Kokkos::View<..., Kokkos::HostSpace> objects and
    ///       opens one FILE*.
    ///
    /// @code
    ///     fierro::Logger log(MPI_COMM_WORLD, fierro::LogLevel::Info, "Fierro_log_");
    /// @endcode
    inline explicit Logger(MPI_Comm    comm           = MPI_COMM_WORLD,
                           LogLevel    level          = LogLevel::Info,
                           const char* file_prefix    = "Fierro_log_",
                           size_t      capacity_bytes = static_cast<size_t>(1) << 20)
        : size_(1)
        , comm_(comm)
        , file_(nullptr)
    {
        int r = 0, sz = 1;
        MPI_Comm_rank(comm_, &r);
        MPI_Comm_size(comm_, &sz);
        size_ = sz;

        // Allocate Kokkos-backed storage. The Views own the memory; the
        // Handle exposes raw pointers into them for kernel capture. The
        // Kokkos runtime must be initialized before the Logger ctor runs
        // (MATAR_INITIALIZE covers this in Fierro's main).
        buf_view_     = Kokkos::View<char*  , Kokkos::HostSpace>("fierro.logger.buf", capacity_bytes);
        used_view_    = Kokkos::View<size_t , Kokkos::HostSpace>("fierro.logger.used");
        dropped_view_ = Kokkos::View<size_t , Kokkos::HostSpace>("fierro.logger.dropped");
        *used_view_.data()    = 0;
        *dropped_view_.data() = 0;

        h_.buf_ptr = buf_view_.data();
        h_.used    = used_view_.data();
        h_.dropped = dropped_view_.data();
        h_.cap     = capacity_bytes;
        h_.rank    = r;
        h_.level   = level;

        if (file_prefix != nullptr) {
            open_file_for_rank_(file_prefix);
        }
    }

    /// @brief Best-effort destructor. If MPI is still initialized it calls
    ///        flush() (collective!) so no messages are lost; then closes the
    ///        file. If MPI has already been finalized, only the rank-local
    ///        buffer is written to the file and the stdout gather is skipped.
    ///
    /// @warning Because ~Logger is collective when MPI is live, every rank
    ///          must destroy its Logger at the same program point. The
    ///          intended pattern is to make Logger a Driver member so its
    ///          destruction coincides with Driver destruction, before
    ///          MPI_Finalize.
    inline ~Logger() {
        int mpi_initialized = 0;
        int mpi_finalized   = 0;
        MPI_Initialized(&mpi_initialized);
        MPI_Finalized(&mpi_finalized);

        if (mpi_initialized && !mpi_finalized) {
            // Collective path: safe to gather to rank 0 stdout.
            flush();
        } else {
            // MPI gone: still write the rank-local buffer to this rank's file.
            drain_local_only_();
        }
        close_file_();
    }

    /// @brief Copying a Logger is forbidden (the owned FILE* and Kokkos::View
    ///        buffers have unique ownership). Use Logger::handle() for cheap
    ///        value-capture into kernels, or pass by reference across solver
    ///        boundaries.
    Logger(const Logger&)            = delete;
    Logger& operator=(const Logger&) = delete;

    // ------------------------------------------------------------------------
    //  Host-side logging convenience (forwards to h_.log / h_.<level>)
    // ------------------------------------------------------------------------

    /// @brief Host-side convenience that forwards to h_.log(lvl, fmt, args...).
    /// @copydetails fierro::Logger::Handle::log
    template<typename... Args>
    inline void log(LogLevel lvl, const char* fmt, Args... args) const {
        h_.log(lvl, fmt, args...);
    }

    /// @brief Host-side convenience for Trace. Non-collective. Host-thread-safe.
    template<typename... Args>
    inline void trace(const char* fmt, Args... args) const { h_.trace(fmt, args...); }

    /// @brief Host-side convenience for Debug. Non-collective. Host-thread-safe.
    template<typename... Args>
    inline void debug(const char* fmt, Args... args) const { h_.debug(fmt, args...); }

    /// @brief Host-side convenience for Info. Non-collective. Host-thread-safe.
    template<typename... Args>
    inline void info(const char* fmt, Args... args) const { h_.info(fmt, args...); }

    /// @brief Host-side convenience for Warn. Non-collective. Host-thread-safe.
    template<typename... Args>
    inline void warn(const char* fmt, Args... args) const { h_.warn(fmt, args...); }

    /// @brief Host-side convenience for Error. Non-collective. Host-thread-safe.
    template<typename... Args>
    inline void error(const char* fmt, Args... args) const { h_.error(fmt, args...); }

    // ------------------------------------------------------------------------
    //  Handle access for kernel capture
    // ------------------------------------------------------------------------

    /// @brief Produce a POD view suitable for capturing by value into a kernel.
    ///
    /// @return A trivially-copyable Handle struct (40 bytes). Cheap: copies
    ///         pointers + counters + scalars.
    ///
    /// @note Host-only.
    /// @note The returned Handle captures the CURRENT values of rank/level.
    ///       If Logger::set_level is called later, existing Handles are
    ///       unaffected (their snapshot of `level` is stale). Refresh with
    ///       log.handle() after changing level.
    ///
    /// @code
    ///     auto lh = log.handle();
    ///     FOR_ALL(i, 0, N, { lh.error("elem=%d bad\n", (int)i); });
    /// @endcode
    inline Handle handle() const noexcept { return h_; }

    // ------------------------------------------------------------------------
    //  File prefix management
    // ------------------------------------------------------------------------

    /// @brief Flush and close the current per-rank log file, then open a new
    ///        one at "<prefix><rank>" (mode "w"). Safe to call at any point
    ///        when no other host thread is logging.
    ///
    /// @param prefix  New filename prefix. Pass nullptr to stop file output
    ///                entirely (the current file is closed and none is
    ///                re-opened).
    ///
    /// @note Non-collective by itself; but the nested flush() IS collective.
    ///       Recommended usage is to call on every rank at the same program
    ///       point (e.g. after parsing YAML in Driver).
    ///
    /// @code
    ///     log.set_file_prefix("run42_");  // -> run42_0, run42_1, ...
    /// @endcode
    inline void set_file_prefix(const char* prefix) {
        // flush current buffer so nothing is lost, then reopen under new name
        flush();
        close_file_();
        if (prefix != nullptr) {
            open_file_for_rank_(prefix);
        }
    }

    // ------------------------------------------------------------------------
    //  The one collective: flush()
    // ------------------------------------------------------------------------

    /// @brief Drain the rank-local buffer. Writes to the rank's own file (no
    ///        MPI) and to rank-0 stdout (via one MPI_Gather + one
    ///        MPI_Gatherv<MPI_CHAR>) in strict rank order. Resets both the
    ///        byte and drop counters to zero.
    ///
    /// @note COLLECTIVE over the Logger's communicator. Must be reached by
    ///       every rank the same number of times. Call at phase boundaries
    ///       (end of initialize/setup/each cycle/execute/finalize). Never
    ///       call from inside a loop whose iteration count varies across
    ///       ranks.
    ///
    /// @note NOT host-thread-safe. The hot-path append IS thread-safe, but
    ///       flush resets the counters concurrently with potential appends.
    ///       Call flush from a single host thread.
    ///
    /// @note If the overflow counter is nonzero, a single "[rank N][WARN ]
    ///       logger: <N> messages dropped (buffer overflow)\n" line is
    ///       emitted on both the file and stdout, then reset.
    ///
    /// @code
    ///     // inside Driver::execute(), at each phase boundary:
    ///     log.flush();
    /// @endcode
    inline void flush() {
        const int rank = h_.rank;
        const int nr   = size_;

        // `used` is allowed to exceed cap (overflow is counted, not stored).
        // Clamp to the actual bytes written.
        const size_t used_now = *h_.used;
        const int    n        = static_cast<int>(used_now < h_.cap ? used_now : h_.cap);

        Kokkos::View<int*, Kokkos::HostSpace> counts, displs;
        if (rank == 0) {
            counts = Kokkos::View<int*, Kokkos::HostSpace>("fierro.logger.counts", nr);
            displs = Kokkos::View<int*, Kokkos::HostSpace>("fierro.logger.displs", nr);
        }
        MPI_Gather(&n, 1, MPI_INT,
                   rank == 0 ? counts.data() : nullptr, 1, MPI_INT,
                   0, comm_);

        int total = 0;
        if (rank == 0) {
            for (int i = 0; i < nr; ++i) {
                displs(i) = total;
                total    += counts(i);
            }
        }

        Kokkos::View<char*, Kokkos::HostSpace> gathered(
            "fierro.logger.gathered",
            rank == 0 ? static_cast<size_t>(total) : static_cast<size_t>(0));

        MPI_Gatherv(h_.buf_ptr, n, MPI_CHAR,
                    rank == 0 ? gathered.data() : nullptr,
                    rank == 0 ? counts.data()   : nullptr,
                    rank == 0 ? displs.data()   : nullptr,
                    MPI_CHAR, 0, comm_);

        if (rank == 0 && total > 0) {
            std::fwrite(gathered.data(), 1, static_cast<size_t>(total), stdout);
            std::fflush(stdout);
        }

        // Every rank: append its own buffer to its own file (no MPI).
        if (file_ != nullptr && n > 0) {
            std::fwrite(h_.buf_ptr, 1, static_cast<size_t>(n), file_);
            std::fflush(file_);
        }

        // Emit a single "dropped" warning line if any messages overflowed.
        const size_t nd = *h_.dropped;
        if (nd > 0) {
            char warn_line[160];
            int w = std::snprintf(
                warn_line, sizeof(warn_line),
                "[rank %d][WARN ] logger: %zu messages dropped (buffer overflow)\n",
                rank, nd);
            if (w > 0) {
                if (file_ != nullptr) {
                    std::fwrite(warn_line, 1, static_cast<size_t>(w), file_);
                    std::fflush(file_);
                }
                if (rank == 0) {
                    std::fwrite(warn_line, 1, static_cast<size_t>(w), stdout);
                    std::fflush(stdout);
                }
            }
            *h_.dropped = 0;
        }

        *h_.used = 0;
    }

    // ------------------------------------------------------------------------
    //  Getters / level control
    // ------------------------------------------------------------------------

    /// @brief MPI rank of this Logger in its communicator. Host-only. noexcept.
    inline int rank() const noexcept { return h_.rank; }

    /// @brief Size of the Logger's communicator. Host-only. noexcept.
    inline int size() const noexcept { return size_; }

    /// @brief Current host-side minimum level.
    inline LogLevel level() const noexcept { return h_.level; }

    /// @brief Change the host-side minimum level. Affects subsequent calls on
    ///        this Logger and on Handles produced AFTER this call.
    /// @warning Does not affect Handles already captured into running kernels.
    inline void set_level(LogLevel l) noexcept { h_.level = l; }

    /// @brief MPI communicator bound at construction.
    inline MPI_Comm comm() const noexcept { return comm_; }

private:
    // ------------------------------------------------------------------------
    //  Private helpers
    // ------------------------------------------------------------------------

    /// @brief Build file_path_ and std::fopen in mode "w" (truncating).
    ///        Stores nullptr in file_ if fopen fails (errors are soft: logging
    ///        continues, file output is simply absent).
    inline void open_file_for_rank_(const char* prefix) {
        // "%s%d" keeps us free of any std::string dependency.
        std::snprintf(file_path_, sizeof(file_path_), "%s%d", prefix, h_.rank);
        file_ = std::fopen(file_path_, "w");
    }

    /// @brief Close the per-rank file if one is open. Idempotent.
    inline void close_file_() {
        if (file_ != nullptr) {
            std::fflush(file_);
            std::fclose(file_);
            file_ = nullptr;
        }
    }

    /// @brief Dtor fallback when MPI has already been finalized. Writes the
    ///        rank-local buffer to the file (and a drop-count line, if any),
    ///        then resets counters. Never touches MPI or stdout.
    inline void drain_local_only_() {
        const size_t used_now = *h_.used;
        const size_t n        = used_now < h_.cap ? used_now : h_.cap;
        if (file_ != nullptr && n > 0) {
            std::fwrite(h_.buf_ptr, 1, n, file_);
            std::fflush(file_);
        }
        const size_t nd = *h_.dropped;
        if (nd > 0 && file_ != nullptr) {
            char warn_line[160];
            int w = std::snprintf(
                warn_line, sizeof(warn_line),
                "[rank %d][WARN ] logger: %zu messages dropped (buffer overflow)\n",
                h_.rank, nd);
            if (w > 0) {
                std::fwrite(warn_line, 1, static_cast<size_t>(w), file_);
                std::fflush(file_);
            }
        }
        *h_.used    = 0;
        *h_.dropped = 0;
    }

    // ------------------------------------------------------------------------
    //  State
    // ------------------------------------------------------------------------

    Handle                                   h_{};             ///< View into owned storage
    Kokkos::View<char*  , Kokkos::HostSpace> buf_view_;        ///< Owns the buffer bytes
    Kokkos::View<size_t , Kokkos::HostSpace> used_view_;       ///< Owns the byte counter
    Kokkos::View<size_t , Kokkos::HostSpace> dropped_view_;    ///< Owns the overflow counter
    int                                      size_ = 1;        ///< MPI comm size
    MPI_Comm                                 comm_ = MPI_COMM_WORLD;
    FILE*                                    file_ = nullptr;  ///< Per-rank log file (nullable)

    /// Fixed-size storage for "<prefix><rank>" (no std::string dependency).
    static constexpr int kMaxPath = 256;
    char                                     file_path_[kMaxPath] = {0};
};

/// @brief Null-safe POD snapshot for Kokkos / FOR_ALL. `nullptr` yields a
///        default Handle for which `if (h)` is false and logging is a no-op.
///        Solvers typically write `const auto log_dev = log_handle();` using
///        Solver::log_handle() instead of spelling this out each time.
inline Logger::Handle log_handle(Logger* p) noexcept {
    return p ? p->handle() : Logger::Handle{};
}

} // namespace fierro


// ============================================================================
//  FLOG_DEV -- optional fused single-printf macro for hot device sites
// ============================================================================

/// @def   FLOG_DEV(handle, LEVEL, fmt, ...)
/// @brief Emit a single Kokkos::printf from inside a FOR_ALL / KOKKOS_LAMBDA.
///        Fuses the "[rank N][LEVEL] " tag into the format string via
///        preprocessor string-literal concatenation, so one printf call
///        produces one whole line (no cross-thread interleaving between tag
///        and body).
///
/// @param handle  A fierro::Logger::Handle captured by value into the kernel.
/// @param LEVEL   Bare identifier: TRACE / DEBUG / INFO / WARN / ERROR.
///                (Preprocessor-stringified; must NOT be a variable.)
/// @param fmt     printf format string literal.
/// @param ...     printf arguments.
///
/// @note Device path only (expands to Kokkos::printf). Safe to call from host
///       code too -- Kokkos::printf falls back to std::printf there.
/// @note Does NOT hit the host buffer. Output appears on stdout via the
///       Kokkos::printf path, but NOT in the per-rank log file.
///
/// @code
///     auto lh = log.handle();
///     FOR_ALL(i, 0, N, {
///         if (elem_gid == probe_gid) {
///             FLOG_DEV(lh, ERROR, "negative density in elem %lu\n", elem_gid);
///         }
///     });
/// @endcode
#define FLOG_DEV(handle, LEVEL, fmt, ...) \
    ::Kokkos::printf("[rank %d][" #LEVEL "] " fmt, (handle).rank, ##__VA_ARGS__)

#endif // FIERRO_LOGGER_H
