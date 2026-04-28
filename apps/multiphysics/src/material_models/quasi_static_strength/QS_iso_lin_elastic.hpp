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
#ifndef QS_ISOTROPIC_LINEAR_ELASTIC_H
#define QS_ISOTROPIC_LINEAR_ELASTIC_H


/////////////////////////////////////////////////////////////////////////////
///
/// \fn HypoPlasticity
///
/// \brief stress model for hypo plasticity response
///
///  This is a material model function for the deviatoric stress tensor
///
/// \param Material pressure
/// \param Material stress
/// \param Global ID for the material
/// \param Material ID for the element
/// \param Material state variables
/// \param Material Sound speed
/// \param Material density
/// \param Material specific internal energy
/// \param Element velocity gradient
/// \param Element nodes IDs in the element
/// \param Node node coordinates
/// \param Noe velocity 
/// \param Element volume
/// \param Time time step
/// \param Time coefficient in the Runge Kutta time integration step
///
/////////////////////////////////////////////////////////////////////////////
namespace QSIsotropicLinearElastic {

    const double fuzz = 1.e-16;
    
    /**
     * @brief Initialize the strength state variables for the material points.
     *
     * @param MaterialPoints_eos_state_vars   State variables for the equation of state at each material point.
     * @param MaterialPoints_strength_state_vars   State variables for the strength model at each material point.
     * @param eos_global_vars   Global variables for the equation of state.
     * @param strength_global_vars   Global variables for the strength model.
     * @param elem_in_mat_elem   Mapping from material points to mesh elements.
     * @param num_material_points   Number of material points for this material.
     * @param mat_id   Material ID.
     */
    static void init_strength_state_vars(
        const DRaggedRightArrayKokkos <double> &MaterialPoints_eos_state_vars,
        const DRaggedRightArrayKokkos <double> &MaterialPoints_strength_state_vars,
        const RaggedRightArrayKokkos <double> &eos_global_vars,
        const RaggedRightArrayKokkos <double> &strength_global_vars,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const size_t num_material_points,
        const size_t mat_id)
    {


    }  // end of init_strength_state_vars


    /**
     * @brief Calculate the deviatoric stress tensor for the hypo-plasticity response.
     *
     * Evolves the deviatoric stress tensor using the Jaumann rate and applies J2 plasticity yield condition.
     *
     * @param GaussPoints_vel_grad   Velocity gradient at Gauss points.
     * @param node_coords   Coordinates of the nodes.
     * @param node_vel   Velocities of the nodes.
     * @param nodes_in_elem   Node indices in the element.
     * @param MaterialPoints_pres   Pressure at material points.
     * @param MaterialPoints_stress   Stress tensor at material points (output).
     * @param MaterialPoints_stress_n0   Stress tensor at material points from previous time step.
     * @param MaterialPoints_sspd   Sound speed at material points.
     * @param MaterialPoints_eos_state_vars   State variables for the equation of state at material points.
     * @param MaterialPoints_strength_state_vars   State variables for the strength model at material points.
     * @param MaterialPoints_den   Density at the material point.
     * @param MaterialPoints_sie   Specific internal energy at the material point.
     * @param MaterialPoints_shear_modulii   Shear modulus at material points.
     * @param elem_in_mat_elem   Mapping from material points to mesh elements.
     * @param eos_global_vars   Global variables for the equation of state.
     * @param strength_global_vars   Global variables for the strength model.
     * @param vol   Volume of the element.
     * @param dt   Time step size.
     * @param rk_alpha   Runge-Kutta time integration coefficient.
     * @param time   Current simulation time.
     * @param cycle   Current simulation cycle.
     * @param MaterialPoints_lid   Local ID of the material point.
     * @param mat_id   Material ID.
     * @param gauss_gid   Global ID of the Gauss point.
     * @param elem_gid   Global ID of the element.
     */
    KOKKOS_FUNCTION
    static void calc_stress(
        const DCArrayKokkos<double>  &GaussPoints_vel_grad,
        const DCArrayKokkos <double> &node_coords,
        const DCArrayKokkos <double> &node_vel,
        const DCArrayKokkos<size_t>  &nodes_in_elem,
        const DRaggedRightArrayKokkos<double>  &MaterialPoints_pres,
        const DRaggedRightArrayKokkos<double>  &MaterialPoints_stress,
        const DRaggedRightArrayKokkos<double>  &MaterialPoints_stress_n0,
        const DRaggedRightArrayKokkos<double>  &MaterialPoints_sspd,
        const DRaggedRightArrayKokkos <double> &MaterialPoints_eos_state_vars,
        const DRaggedRightArrayKokkos <double> &MaterialPoints_strength_state_vars,
        const double MaterialPoints_den,
        const double MaterialPoints_sie,
        const DRaggedRightArrayKokkos<double>& MaterialPoints_shear_modulii,
        const DRaggedRightArrayKokkos<size_t>& elem_in_mat_elem,
        const RaggedRightArrayKokkos <double> &eos_global_vars,
        const RaggedRightArrayKokkos <double> &strength_global_vars,
        const double vol,
        const double dt,
        const double rk_alpha,
        const double time,
        const size_t cycle,
        const size_t MaterialPoints_lid,
        const size_t mat_id,
        const size_t gauss_gid,
        const size_t elem_gid)
    {


    } // end of user mat

    KOKKOS_FUNCTION
    static void fill_C_matrix(const RaggedRightArrayKokkos <double> &strength_global_vars,
        double C[6][6],
        const size_t mat_id)
    {
        // setting youngs modulus and poissons ratio
        double E = strength_global_vars(mat_id,0);
        double v = strength_global_vars(mat_id,1);

        // leading coefficient
        double coef = E / ((1-2*v)*(1+v));

        // filling C
        C[0][0] = coef*(1-v);
        C[0][1] = coef*v;
        C[0][2] = coef*v;
        C[0][3] = 0;
        C[0][4] = 0;
        C[0][5] = 0;

        C[1][0] = coef*v;
        C[1][1] = coef*(1-v);
        C[1][2] = coef*v;
        C[1][3] = 0;
        C[1][4] = 0;
        C[1][5] = 0;

        C[2][0] = coef*v;
        C[2][1] = coef*v;
        C[2][2] = coef*(1-v);
        C[2][3] = 0;
        C[2][4] = 0;
        C[2][5] = 0;

        C[3][0] = 0;
        C[3][1] = 0;
        C[3][2] = 0;
        C[3][3] = coef*(1-2*v);
        C[3][4] = 0;
        C[3][5] = 0;

        C[4][0] = 0;
        C[4][1] = 0;
        C[4][2] = 0;
        C[4][3] = 0;
        C[4][4] = coef*(1-2*v);
        C[4][5] = 0;

        C[5][0] = 0;
        C[5][1] = 0;
        C[5][2] = 0;
        C[5][3] = 0;
        C[5][4] = 0;
        C[5][5] = coef*(1-2*v);

    } // end of fill_C_matrix

} // end namespace


#endif // end Header Guard