/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_C_H
#define HEFFTE_C_H

#ifdef __cplusplus
#error "do not use heffte_c.h with a C++ compiler"
#endif
/*
 * C interface for HeFFTe
 */
#include <mpi.h>
#include "heffte_c_defines.h"

/*!
 * \ingroup cfft
 * \brief Populates the fields of the options with the ones default for the backend.
 *
 * \returns Heffte_SUCCESS on success, 1 on invalid backend.
 */
int heffte_set_default_options(int backend, heffte_plan_options *options);

/*!
 * \ingroup cfft
 * \brief Creates a new heffte_plan with the given backend, boxes, communicator, and options.
 *
 * Creates a new heffte_plan with the given set of parameters.
 * \param backend must be Heffte_BACKEND_FFTW, Heffte_BACKEND_MKL or Heffte_BACKEND_CUFFT
 *                and the corresponding backend must be enabled in CMake
 * \param inbox_low is the low indexes of the inbox in heffte::fft3d()
 * \param inbox_high is the high indexes of the inbox in heffte::fft3d()
 * \param inbox_order is the order indexes of the inbox in heffte::fft3d(),
 *                    if set to NULL the default order (0, 1, 2) will be used
 * \param outbox_low is the low indexes of the outbox in heffte::fft3d()
 * \param outbox_high is the high indexes of the outbox in heffte::fft3d()
 * \param outbox_order is the order indexes of the outbox in heffte::fft3d(),
 *                     if set to NULL the default order (0, 1, 2) will be used
 * \param comm is the MPI communicator holding all the boxes
 * \param options is either NULL, which will force the use of the default backend options,
 *                of an instance of heffte_plan_options that will translate
 *                to the options used by the heffte::fft3d()
 * \param plan will be overwritten with a new heffte_plan
 *             must be destroyed by heffte_plan_destroy()
 *
 * \returns Heffte_SUCCESS on success, 1 on invalid backend, 2 on internal exception,
 *          e.g., due to incorrect box configuration, see heffte::fft3d
 */
int heffte_plan_create(int backend, int const inbox_low[3], int const inbox_high[3], int const *inbox_order,
                       int const outbox_low[3], int const outbox_high[3], int const *outbox_order,
                       MPI_Comm const comm, heffte_plan_options const *options, heffte_plan *plan);

/*!
 * \ingroup cfft
 * \brief Creates a new heffte_plan for an r2c operation with the given backend, boxes, communicator, and options.
 *
 * Identical to heffte_plan_create() with the addition of the r2c_direction
 * and the internal class will be heffte::fft3d_r2c
 */
int heffte_plan_create_r2c(int backend, int const inbox_low[3], int const inbox_high[3], int const *inbox_order,
                           int const outbox_low[3], int const outbox_high[3], int const *outbox_order,
                           int r2c_direction, MPI_Comm const comm, heffte_plan_options const *options, heffte_plan *plan);

/*!
 * \ingroup cfft
 * \brief Destory a heffte_plan, call the destructor of the internal object.
 *
 * \param plan to be destoryed
 *
 * \returns Heffte_SUCCESS on success, 3 on corrupted object
 */
int heffte_plan_destroy(heffte_plan plan);

/*!
 * \ingroup cfft
 * \brief Wrapper around heffte::fft3d::size_inbox() and heffte::fft3d_r2c::size_inbox()
 */
int heffte_size_inbox(heffte_plan const plan);
/*!
 * \ingroup cfft
 * \brief Wrapper around heffte::fft3d::size_outbox() and heffte::fft3d_r2c::size_outbox()
 */
int heffte_size_outbox(heffte_plan const plan);
/*!
 * \ingroup cfft
 * \brief Wrapper around heffte::fft3d::size_workspace() and heffte::fft3d_r2c::size_workspace()
 */
int heffte_size_workspace(heffte_plan const plan);

/*!
 * \ingroup cfft
 * \brief Returns the backend used in the create command.
 */
int heffte_get_backend(heffte_plan const plan);
/*!
 * \ingroup cfft
 * \brief Returns 1 if heffte_plan_create_r2c() was called and 0 otherwise.
 */
int heffte_is_r2c(heffte_plan const plan);

/*!
 * \ingroup cfft
 * \brief Wrapper around heffte::fft3d::forward() and heffte::fft3d_r2c::forward()
 *
 * Performs forward Fourier transform from real input data to complex output data,
 * both in single precision.
 * \param plan is an already created heffte_plan
 * \param input has size at least heffte_size_inbox() and contains the input data
 * \param output has size at least heffte_size_outbox() of complex-single precision numbers
 * \param scale is Heffte_SCALE_NONE, Heffte_SCALE_FULL or Heffte_SCALE_SYMMETRIC,
 *              see heffte::scale for details
 *
 * The \b _s2c suffix indicates the real/complex types and precision using standard BLAS
 * conventions:
 * <table>
 * <tr><td> Input </td><td> Output </td><td> Method </td></tr>
 * <tr><td> float </td><td> std::complex<float> </td><td> heffte_forward_s2c() </td></tr>
 * <tr><td> std::complex<float> </td><td> std::complex<float> </td><td> heffte_forward_c2c() </td></tr>
 * <tr><td> double </td><td> std::complex<double> </td><td> heffte_forward_d2z() </td></tr>
 * <tr><td> std::complex<double> </td><td> std::complex<double> </td><td> heffte_forward_z2z() </td></tr>
 * </table>
 */
void heffte_forward_s2c(heffte_plan const plan, float const *input, void *output, int scale);
/*!
 * \ingroup cfft
 * \brief Forward transform, single precision, complex to complex, see heffte_forward_s2c()
 */
void heffte_forward_c2c(heffte_plan const plan, void const *input, void *output, int scale);
/*!
 * \ingroup cfft
 * \brief Forward transform, double precision, real to complex, see heffte_forward_s2c()
 */
void heffte_forward_d2z(heffte_plan const plan, double const *input, void *output, int scale);
/*!
 * \ingroup cfft
 * \brief Forward transform, double precision, complex to complex, see heffte_forward_s2c()
 */
void heffte_forward_z2z(heffte_plan const plan, void const *input, void *output, int scale);

/*!
 * \ingroup cfft
 * \brief Wrapper around heffte::fft3d::forward() and heffte::fft3d_r2c::forward() with user allocated workspace.
 *
 * Accepts a user allocated workspace buffer, but otherwise identical to heffte_forward_s2c()
 * \param plan is an already created heffte_plan
 * \param input has size at least heffte_size_inbox() and contains the input data
 * \param output has size at least heffte_size_outbox() of complex-single precision numbers
 * \param workspace has size at least heffte_size_workspace() of complex-single precision numbers
 * \param scale is Heffte_SCALE_NONE, Heffte_SCALE_FULL or Heffte_SCALE_SYMMETRIC,
 *              see heffte::scale for details
 *
 * The method has four variants corresponding to different precision and real/complex types:
 * <table>
 * <tr><td> Input </td><td> Output </td><td> Workspace </td><td> Method </td></tr>
 * <tr><td> float </td><td> std::complex<float> </td><td> std::complex<float> </td><td> heffte_forward_s2c_buffered() </td></tr>
 * <tr><td> std::complex<float> </td><td> std::complex<float> </td><td> std::complex<float> </td><td> heffte_forward_c2c_buffered() </td></tr>
 * <tr><td> double </td><td> std::complex<double> </td><td> std::complex<double> </td><td> heffte_forward_d2z_buffered() </td></tr>
 * <tr><td> std::complex<double> </td><td> std::complex<double> </td><td> std::complex<double> </td><td> heffte_forward_z2z_buffered() </td></tr>
 * </table>
 */
void heffte_forward_s2c_buffered(heffte_plan const plan, float const *input, void *output, void *workspace, int scale);
/*!
 * \ingroup cfft
 * \brief Forward transform, buffered, double precision, complex to complex, see heffte_forward_s2c_buffered()
 */
void heffte_forward_c2c_buffered(heffte_plan const plan, void const *input, void *output, void *workspace, int scale);
/*!
 * \ingroup cfft
 * \brief Forward transform, buffered, double precision, complex to complex, see heffte_forward_s2c_buffered()
 */
void heffte_forward_d2z_buffered(heffte_plan const plan, double const *input, void *output, void *workspace, int scale);
/*!
 * \ingroup cfft
 * \brief Forward transform, buffered, double precision, complex to complex, see heffte_forward_s2c_buffered()
 */
void heffte_forward_z2z_buffered(heffte_plan const plan, void const *input, void *output, void *workspace, int scale);

/*!
 * \ingroup cfft
 * \brief Wrapper around heffte::fft3d::backward() and heffte::fft3d_r2c::backward()
 *
 * Performs backward Fourier transform from complex input data to real output data,
 * both in single precision.
 * \param plan is an already created heffte_plan
 * \param input has size at least heffte_size_inbox() of complex-single precision numbers
 * \param output has size at least heffte_size_outbox() will be overwritten with the output
 * \param scale is Heffte_SCALE_NONE, Heffte_SCALE_FULL or Heffte_SCALE_SYMMETRIC,
 *              see heffte::scale for details
 *
 * The \b _s2c suffix indicates the real/complex types and precision using standard BLAS
 * conventions:
 * <table>
 * <tr><td> Input </td><td> Output </td><td> Method </td></tr>
 * <tr><td> std::complex<float> </td><td> float </td><td> heffte_backward_c2s() </td></tr>
 * <tr><td> std::complex<float> </td><td> std::complex<float> </td><td> heffte_backward_c2c() </td></tr>
 * <tr><td> std::complex<double> </td><td> double </td><td> heffte_backward_z2d() </td></tr>
 * <tr><td> std::complex<double> </td><td> std::complex<double> </td><td> heffte_backward_z2z() </td></tr>
 * </table>
 */
void heffte_backward_c2s(heffte_plan const plan, void const *input, float *output, int scale);
/*!
 * \ingroup cfft
 * \brief Backward transform, single precision, complex to complex, see heffte_backward_c2s()
 */
void heffte_backward_c2c(heffte_plan const plan, void const *input, void *output, int scale);
/*!
 * \ingroup cfft
 * \brief Backward transform, double precision, complex to real, see heffte_backward_c2s()
 */
void heffte_backward_z2d(heffte_plan const plan, void const *input, double *output, int scale);
/*!
 * \ingroup cfft
 * \brief Backward transform, double precision, complex to complex, see heffte_backward_c2s()
 */
void heffte_backward_z2z(heffte_plan const plan, void const *input, void *output, int scale);

/*!
 * \ingroup cfft
 * \brief Wrapper around heffte::fft3d::backward() and heffte::fft3d_r2c::backward() with user allocated workspace.
 *
 * Accepts a user allocated workspace buffer, but otherwise identical to heffte_backward_s2c()
 * \param plan is an already created heffte_plan
 * \param input has size at least heffte_size_inbox() of complex-single precision numbers
 * \param output has size at least heffte_size_outbox() will be overwritten with the output
 * \param workspace has size at least heffte_size_workspace() of complex-single precision numbers
 * \param scale is Heffte_SCALE_NONE, Heffte_SCALE_FULL or Heffte_SCALE_SYMMETRIC,
 *              see heffte::scale for details
 *
 * The method has four variants corresponding to different precision and real/complex types:
 * <table>
 * <tr><td> Input </td><td> Output </td><td> Workspace </td><td> Method </td></tr>
 * <tr><td> std::complex<float> </td><td> float </td><td> std::complex<float> </td><td> heffte_backward_c2s_buffered() </td></tr>
 * <tr><td> std::complex<float> </td><td> std::complex<float> </td><td> std::complex<float> </td><td> heffte_backward_c2c_buffered() </td></tr>
 * <tr><td> std::complex<double> </td><td> double </td><td> std::complex<double> </td><td> heffte_backward_z2d_buffered() </td></tr>
 * <tr><td> std::complex<double> </td><td> std::complex<double> </td><td> std::complex<double> </td><td> heffte_backward_z2z_buffered() </td></tr>
 * </table>
 */
void heffte_backward_c2s_buffered(heffte_plan const plan, void const *input, float *output, void *workspace, int scale);
/*!
 * \ingroup cfft
 * \brief Backward transform, buffered, single precision, complex to complex, see heffte_backward_c2s_buffered()
 */
void heffte_backward_c2c_buffered(heffte_plan const plan, void const *input, void *output, void *workspace, int scale);
/*!
 * \ingroup cfft
 * \brief Backward transform, buffered, double precision, complex to real, see heffte_backward_c2s_buffered()
 */
void heffte_backward_z2d_buffered(heffte_plan const plan, void const *input, double *output, void *workspace, int scale);
/*!
 * \ingroup cfft
 * \brief Backward transform, buffered, double precision, complex to complex, see heffte_backward_c2s_buffered()
 */
void heffte_backward_z2z_buffered(heffte_plan const plan, void const *input, void *output, void *workspace, int scale);

#endif /* HEFFTE_C_H */
