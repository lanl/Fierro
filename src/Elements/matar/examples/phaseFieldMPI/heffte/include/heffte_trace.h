/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_TRACE_H
#define HEFFTE_TRACE_H

#include "heffte_utils.h"

/*!
 * \ingroup fft3d
 * \addtogroup hefftetrace Event tracing
 *
 * HeFFTe provides a simple tracing capability that will log different events
 * in both a forward and backward transform.
 * Tracing is disabled by default (in order to remove any overhead),
 * and is designed to have minimal overhead at run time.
 *
 * \anchor hefftetracedescription
 * Tracing must first be enabled in CMake with
 * \code
 *  -D Heffte_ENABLE_TRACING=ON
 * \endcode
 * and in the client code must be initialized with
 * \code
 *  init_tracing("root-file-name");
 * \endcode
 * which must be called after MPI_Init().
 * Then, any FFT transform will log different events, e.g., packing, unpacking,
 * communicating all-to-all, etc. Each event will be times beginning to end
 * using MPI_Wtime().
 * Tracing is finalized with
 * \code
 *  finalize_tracing();
 * \endcode
 * which must be called before MPI_Finalize().
 * The finalize call will write a number of text files (one per MPI rank)
 * each with name "root-file-name_rank_number.txt"
 * that will hold the list of events encountered during the calls to HeFFTe.
 *
 * \b Note: if tracing is not enabled, all tracing methods are no-op.
 */

namespace heffte {

    #ifdef Heffte_ENABLE_TRACING
    /*!
     * \ingroup hefftetrace
     * \brief A tracing event.
     *
     * The events that are being traced will be associated with start time,
     * duration, and name of the event.
     */
    struct event {
        //! \brief Name of the event.
        std::string name;
        //! \brief Start time according to MPI high-precision clock.
        double start_time;
        //! \brief Duration according to the MPI high-precision clock.
        double duration;
    };

    /*!
     * \ingroup hefftetrace
     * \brief Logs the list of events (declared in heffte_reshape3d.cpp).
     */
    extern std::deque<event> event_log;
    /*!
     * \ingroup hefftetrace
     * \brief Root filename to write out the traces (decaled in heffte_reshape3d.cpp).
     */
    extern std::string log_filename;

    /*!
     * \brief hefftetrace
     * \brief Declares an event to be traced.
     *
     * The constructor will record the current wall-clock-time and the name of the event,
     * the destructor will measure the wall-clock-time and will add the event to the event_log
     * saving the name, start time and duration.
     */
    struct add_trace {
        //! \brief Defines a trace.
        add_trace(std::string s) : name(s), start_time(MPI_Wtime()){}
        //! \brief Saves a trace.
        ~add_trace(){
            double duration = MPI_Wtime() - start_time;
            event_log.push_back({std::move(name), start_time, duration});
        }
        //! \brief Defines the name of a trace.
        std::string name;
        //! \brief Saves the starting wall-clock-time of the trace.
        double start_time;
    };

    /*!
     * \ingroup hefftetrace
     * \brief Initialize tracing and remember the root filename for output, see the \ref hefftetracedescription "Detailed Description"
     */
    inline void init_tracing(std::string root_filename){
        event_log = std::deque<event>();
        log_filename = root_filename + "_" + std::to_string(mpi::world_rank()) + ".log";
    }

    /*!
     * \ingroup hefftetrace
     * \brief Finalize tracing and write the result to a file, see the \ref hefftetracedescription "Detailed Description"
     */
    inline void finalize_tracing(){
        std::ofstream ofs(log_filename);
        ofs.precision(12);

        for(auto const &e : event_log)
            ofs << std::setw(40) << e.name << std::setw(20) << e.start_time << std::setw(20) << e.duration << "\n";
    }

    #else

    struct add_trace{
        add_trace(std::string){}
    };
    inline void init_tracing(std::string){}
    inline void finalize_tracing(){}

    #endif

}

#endif        //  #ifndef HEFFTE_TRACE_H
