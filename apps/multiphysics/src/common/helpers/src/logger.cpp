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

#include "logger.hpp"

namespace fierro {

Logger::Logger(MPI_Comm    comm,
                           LogLevel    level,
                           const char* file_prefix,
                           size_t      capacity_bytes) : size_(1) , comm_(comm) , file_(nullptr)
{
    int r = 0, sz = 1;
    MPI_Comm_rank(comm_, &r);
    MPI_Comm_size(comm_, &sz);
    size_ = sz;

    const size_t cap_eff    = (capacity_bytes < 2u) ? 2u : capacity_bytes;
    const size_t file_cap   = cap_eff / 2;
    const size_t stdout_cap = cap_eff - file_cap;

    buf_view_ = Kokkos::View<char*, Kokkos::HostSpace>("fierro.logger.buf", cap_eff);
    used_view_         = Kokkos::View<size_t, Kokkos::HostSpace>("fierro.logger.used");
    dropped_view_      = Kokkos::View<size_t, Kokkos::HostSpace>("fierro.logger.dropped");
    stdout_used_view_  = Kokkos::View<size_t, Kokkos::HostSpace>("fierro.logger.stdout_used");
    stdout_dropped_view_ =
        Kokkos::View<size_t, Kokkos::HostSpace>("fierro.logger.stdout_dropped");
    *used_view_.data()           = 0;
    *dropped_view_.data()         = 0;
    *stdout_used_view_.data()    = 0;
    *stdout_dropped_view_.data() = 0;

    char* const base = buf_view_.data();
    h_.buf_ptr         = base;
    h_.used            = used_view_.data();
    h_.dropped         = dropped_view_.data();
    h_.cap             = file_cap;
    h_.stdout_buf_ptr  = base + file_cap;
    h_.stdout_used     = stdout_used_view_.data();
    h_.stdout_dropped  = stdout_dropped_view_.data();
    h_.stdout_cap      = stdout_cap;
    h_.rank            = r;
    h_.level           = level;

    if (file_prefix != nullptr) {
        open_file_for_rank_(file_prefix);
    }
}

Logger::~Logger()
{
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

void Logger::set_file_prefix(const char* prefix)
{
    // flush current buffer so nothing is lost, then reopen under new name
    flush();
    close_file_();
    if (prefix != nullptr) {
        open_file_for_rank_(prefix);
    }
}

void Logger::flush()
{
    const int rank = h_.rank;
    const int nr   = size_;

    // stdout ring -> rank 0 terminal (MPI gather)
    {
        const size_t used_stdout = *h_.stdout_used;
        const int    n_out       = static_cast<int>(
            used_stdout < h_.stdout_cap ? used_stdout : h_.stdout_cap);

        Kokkos::View<int*, Kokkos::HostSpace> counts, displs;
        if (rank == 0) {
            counts = Kokkos::View<int*, Kokkos::HostSpace>("fierro.logger.sout_counts", nr);
            displs = Kokkos::View<int*, Kokkos::HostSpace>("fierro.logger.sout_displs", nr);
        }
        MPI_Gather(&n_out, 1, MPI_INT,
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
            "fierro.logger.stdout_gathered",
            rank == 0 ? static_cast<size_t>(total) : static_cast<size_t>(0));

        MPI_Gatherv(h_.stdout_buf_ptr, n_out, MPI_CHAR,
                    rank == 0 ? gathered.data() : nullptr,
                    rank == 0 ? counts.data()   : nullptr,
                    rank == 0 ? displs.data()   : nullptr,
                    MPI_CHAR, 0, comm_);

        if (rank == 0 && total > 0) {
            std::fwrite(gathered.data(), 1, static_cast<size_t>(total), stdout);
            std::fflush(stdout);
        }
    }

    // file ring -> this rank's log file (local only)
    {
        const size_t used_now = *h_.used;
        const int    n      = static_cast<int>(used_now < h_.cap ? used_now : h_.cap);
        if (file_ != nullptr && n > 0) {
            std::fwrite(h_.buf_ptr, 1, static_cast<size_t>(n), file_);
            std::fflush(file_);
        }

        const size_t nd = *h_.dropped;
        if (nd > 0) {
            char warn_line[200];
            int w = std::snprintf(
                warn_line, sizeof(warn_line),
                "[rank %d][WARN ] logger (file pool): %zu messages dropped "
                "(buffer overflow)\n",
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

    const size_t nsd = *h_.stdout_dropped;
    if (nsd > 0) {
        char warn_line[200];
        int w = std::snprintf(
            warn_line, sizeof(warn_line),
            "[rank %d][WARN ] logger (stdout pool): %zu messages dropped "
            "(buffer overflow)\n",
            rank, nsd);
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
        *h_.stdout_dropped = 0;
    }
    *h_.stdout_used = 0;
}

void Logger::open_file_for_rank_(const char* prefix)
{
    // "%s%d" keeps us free of any std::string dependency.
    std::snprintf(file_path_, sizeof(file_path_), "%s%d", prefix, h_.rank);
    file_ = std::fopen(file_path_, "w");
}

void Logger::close_file_()
{
    if (file_ != nullptr) {
        std::fflush(file_);
        std::fclose(file_);
        file_ = nullptr;
    }
}

void Logger::drain_local_only_()
{
    const size_t used_file = *h_.used;
    const size_t n_file    = used_file < h_.cap ? used_file : h_.cap;
    if (file_ != nullptr && n_file > 0) {
        std::fwrite(h_.buf_ptr, 1, n_file, file_);
        std::fflush(file_);
    }
    const size_t used_out = *h_.stdout_used;
    const size_t n_out    = used_out < h_.stdout_cap ? used_out : h_.stdout_cap;
    if (file_ != nullptr && n_out > 0 && h_.stdout_buf_ptr != nullptr) {
        std::fwrite(h_.stdout_buf_ptr, 1, n_out, file_);
        std::fflush(file_);
    }
    const size_t nd = *h_.dropped;
    if (nd > 0 && file_ != nullptr) {
        char warn_line[200];
        int w = std::snprintf(
            warn_line, sizeof(warn_line),
            "[rank %d][WARN ] logger (file pool): %zu messages dropped "
            "(buffer overflow)\n",
            h_.rank, nd);
        if (w > 0) {
            std::fwrite(warn_line, 1, static_cast<size_t>(w), file_);
            std::fflush(file_);
        }
    }
    const size_t nsd = *h_.stdout_dropped;
    if (nsd > 0 && file_ != nullptr) {
        char warn_line[200];
        int w = std::snprintf(
            warn_line, sizeof(warn_line),
            "[rank %d][WARN ] logger (stdout pool): %zu messages dropped "
            "(buffer overflow)\n",
            h_.rank, nsd);
        if (w > 0) {
            std::fwrite(warn_line, 1, static_cast<size_t>(w), file_);
            std::fflush(file_);
        }
    }
    *h_.used            = 0;
    *h_.dropped         = 0;
    *h_.stdout_used     = 0;
    *h_.stdout_dropped  = 0;
}

} // namespace fierro
