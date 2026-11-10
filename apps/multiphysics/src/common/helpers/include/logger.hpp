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
/// @brief  Simple MPI-aware, host-only printf-style logger for Fierro.
///
/// Every rank writes its lines to its own file "<file_prefix><rank>"
/// (default Fierro_log_0, Fierro_log_1, ...). Rank 0 also writes to stdout.
/// Lines are tagged "[rank N][LEVEL] ".
///
///     log.info("num_elems=%zu\n", n);
///     log.set_level(fierro::LogLevel::Warn);   // drop trace/debug/info
///
/// No method is collective: logging and flush() never call MPI, so they are
/// safe inside loops whose trip count varies across ranks. Host-only; do not
/// call from inside FOR_ALL / KOKKOS_LAMBDA.

#ifndef FIERRO_LOGGER_H
#define FIERRO_LOGGER_H

#include <mpi.h>
#include <cstdio>
#include <cstdarg>

#if defined(__GNUC__) || defined(__clang__)
#  define FIERRO_PRINTF_LIKE(fmt_idx, args_idx) \
       __attribute__((format(printf, fmt_idx, args_idx)))
#else
#  define FIERRO_PRINTF_LIKE(fmt_idx, args_idx)
#endif

namespace fierro {

/// @brief Severity ordering: Trace < Debug < Info < Warn < Error < Off.
///        Messages below Logger::level() are dropped; Off disables all output.
enum class LogLevel : int {
    Trace = 0,
    Debug = 1,
    Info  = 2,
    Warn  = 3,
    Error = 4,
    Off   = 5
};

class Logger {
public:
    /// @param comm         Communicator used to determine rank/size. MPI must
    ///                     already be initialized.
    /// @param level        Initial minimum level.
    /// @param file_prefix  Per-rank file prefix; nullptr disables file output.
    explicit Logger(MPI_Comm    comm        = MPI_COMM_WORLD,
                    LogLevel    level       = LogLevel::Info,
                    const char* file_prefix = "Fierro_log_");
    ~Logger();

    Logger(const Logger&)            = delete;
    Logger& operator=(const Logger&) = delete;

    void trace(const char* fmt, ...) const FIERRO_PRINTF_LIKE(2, 3);
    void debug(const char* fmt, ...) const FIERRO_PRINTF_LIKE(2, 3);
    void info (const char* fmt, ...) const FIERRO_PRINTF_LIKE(2, 3);
    void warn (const char* fmt, ...) const FIERRO_PRINTF_LIKE(2, 3);
    void error(const char* fmt, ...) const FIERRO_PRINTF_LIKE(2, 3);
    void log(LogLevel lvl, const char* fmt, ...) const FIERRO_PRINTF_LIKE(3, 4);

    /// @brief Flush this rank's file (and stdout on rank 0). Not collective.
    void flush() const;

    int      rank()  const noexcept { return rank_; }
    int      size()  const noexcept { return size_; }
    LogLevel level() const noexcept { return level_; }
    void     set_level(LogLevel l) noexcept { level_ = l; }
    MPI_Comm comm()  const noexcept { return comm_; }

private:
    void vlog_(LogLevel lvl, const char* fmt, va_list ap) const;

    int      rank_  = 0;
    int      size_  = 1;
    LogLevel level_ = LogLevel::Info;
    MPI_Comm comm_  = MPI_COMM_WORLD;
    FILE*    file_  = nullptr;
};

} // namespace fierro

#endif // FIERRO_LOGGER_H
