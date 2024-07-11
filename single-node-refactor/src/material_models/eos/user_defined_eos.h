/**********************************************************************************************
� 2020. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and
to permit others to do so.
This program is open source under the BSD-3 License.
Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:
1.  Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer.
2.  Redistributions in binary form must reproduce the above copyright notice, this list of
conditions and the following disclaimer in the documentation and/or other materials
provided with the distribution.
3.  Neither the name of the copyright holder nor the names of its contributors may be used
to endorse or promote products derived from this software without specific prior
written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**********************************************************************************************/

#ifndef USER_DEFINED_EOS_H
#define USER_DEFINED_EOS_H

/////////////////////////////////////////////////////////////////////////////
///
/// \fn UserDefinedEOSModel
///
/// \brief user defined EOS model
///
/// This is the user material model function for the equation of state
/// An eos function must be supplied or the code will fail to run.
/// The pressure and sound speed can be calculated from an analytic eos.
/// The pressure can also be calculated using p = -1/3 Trace(Stress)
///
/// \param Element pressure
/// \param Element stress
/// \param Global ID for the element
/// \param Material ID for the element
/// \param Element state variables
/// \param Element Sound speed
/// \param Material density
/// \param Material specific internal energy
///
/////////////////////////////////////////////////////////////////////////////
namespace UserDefinedEOSModel
{

    KOKKOS_FUNCTION
    static void calc_pressure(const DCArrayKokkos<double>& elem_pres,
                              const DCArrayKokkos<double>& elem_stress,
                              const size_t elem_gid,
                              const size_t mat_id,
                              const DCArrayKokkos<double>& elem_state_vars,
                              const DCArrayKokkos<double>& elem_sspd,
                              const double den,
                              const double sie,
                              const RaggedRightArrayKokkos<double> &eos_global_vars)
    {
        // -----------------------------------------------------------------------------
        // Required variables are here
        // ------------------------------------------------------------------------------

        // -----------------------------------------------------------------------------
        // The user must coding goes here
        // ------------------------------------------------------------------------------

        return;
    } // end for user_eos_model

    KOKKOS_FUNCTION
    static void calc_sound_speed(const DCArrayKokkos<double>& elem_pres,
                                 const DCArrayKokkos<double>& elem_stress,
                                 const size_t elem_gid,
                                 const size_t mat_id,
                                 const DCArrayKokkos<double>& elem_state_vars,
                                 const DCArrayKokkos<double>& elem_sspd,
                                 const double den,
                                 const double sie,
                                 const RaggedRightArrayKokkos<double> &eos_global_vars)
    {
        
        // -----------------------------------------------------------------------------
        // Required variables are here
        // ------------------------------------------------------------------------------

        // -----------------------------------------------------------------------------
        // The user must coding goes here
        // ------------------------------------------------------------------------------

        return;
    } // end func

} // end namespace



/////////////////////////////////////////////////////////////////////////////
///
/// \fn fcn_name
///
/// \brief <insert brief description>
///
/// <Insert longer more detailed description which
/// can span multiple lines if needed>
///
/// \param <function parameter description>
/// \param <function parameter description>
/// \param <function parameter description>
///
/// \return <return type and definition description if not void>
///
/////////////////////////////////////////////////////////////////////////////
// -----------------------------------------------------------------------------
// This is a place holder EOS interface to add another user defined EOS
// ------------------------------------------------------------------------------
namespace NotionalEOSModel {

    KOKKOS_FUNCTION
    static void calc_pressure(const DCArrayKokkos<double>& elem_pres,
        const DCArrayKokkos<double>& elem_stress,
        const size_t elem_gid,
        const size_t mat_id,
        const DCArrayKokkos<double>& elem_state_vars,
        const DCArrayKokkos<double>& elem_sspd,
        const double den,
        const double sie,
        const RaggedRightArrayKokkos<double> &eos_global_vars)
    {
        // pressure of a void is 0
        elem_pres(elem_gid) = 0.0;

        return;
    } // end func

    KOKKOS_FUNCTION
    static void calc_sound_speed(const DCArrayKokkos<double>& elem_pres,
        const DCArrayKokkos<double>& elem_stress,
        const size_t elem_gid,
        const size_t mat_id,
        const DCArrayKokkos<double>& elem_state_vars,
        const DCArrayKokkos<double>& elem_sspd,
        const double den,
        const double sie,
        const RaggedRightArrayKokkos<double> &eos_global_vars)
    {

        // sound speed of a void is 0, machine small must be used for CFL calculation
        elem_sspd(elem_gid) = 1.0e-32;

        return;
    } // end func

} // end namespace


#endif // end Header Guard