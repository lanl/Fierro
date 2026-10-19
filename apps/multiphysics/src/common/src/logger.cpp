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

namespace {
const char* level_tag(LogLevel l)
{
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
} // namespace

Logger::Logger(MPI_Comm comm, LogLevel level, const char* file_prefix)
    : level_(level)
    , comm_(comm)
{
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);

    if (file_prefix != nullptr) {
        char path[256];
        std::snprintf(path, sizeof(path), "%s%d", file_prefix, rank_);
        file_ = std::fopen(path, "w");
    }
}

Logger::~Logger()
{
    flush();
    if (file_ != nullptr) {
        std::fclose(file_);
    }
}

void Logger::vlog_(LogLevel lvl, const char* fmt, va_list ap) const
{
    if (static_cast<int>(lvl) < static_cast<int>(level_)) return;

    if (file_ != nullptr) {
        va_list ap_file;
        va_copy(ap_file, ap);
        std::fprintf(file_, "[rank %d][%s] ", rank_, level_tag(lvl));
        std::vfprintf(file_, fmt, ap_file);
        va_end(ap_file);
    }
    if (rank_ == 0) {
        va_list ap_out;
        va_copy(ap_out, ap);
        std::fprintf(stdout, "[rank %d][%s] ", rank_, level_tag(lvl));
        std::vfprintf(stdout, fmt, ap_out);
        va_end(ap_out);
    }
}

void Logger::log(LogLevel lvl, const char* fmt, ...) const
{
    va_list ap;
    va_start(ap, fmt);
    vlog_(lvl, fmt, ap);
    va_end(ap);
}

#define FIERRO_LOGGER_LEVEL_METHOD(name, lvl)         \
    void Logger::name(const char* fmt, ...) const     \
    {                                                 \
        va_list ap;                                   \
        va_start(ap, fmt);                            \
        vlog_(lvl, fmt, ap);                          \
        va_end(ap);                                   \
    }

FIERRO_LOGGER_LEVEL_METHOD(trace, LogLevel::Trace)
FIERRO_LOGGER_LEVEL_METHOD(debug, LogLevel::Debug)
FIERRO_LOGGER_LEVEL_METHOD(info,  LogLevel::Info)
FIERRO_LOGGER_LEVEL_METHOD(warn,  LogLevel::Warn)
FIERRO_LOGGER_LEVEL_METHOD(error, LogLevel::Error)

#undef FIERRO_LOGGER_LEVEL_METHOD

void Logger::flush() const
{
    if (file_ != nullptr) std::fflush(file_);
    if (rank_ == 0) std::fflush(stdout);
}

} // namespace fierro
