#ifndef RDH_H
#define RDH_H

#include "matar.h"
#include "mesh.h"
#include "state.h"
#include "ref_elem.h"
#include "ref_surf_elem.h"
#include <cmath>

#ifndef PI
#define PI 3.14159265358979323846
#endif

using namespace mtr;

void read_mesh_ensight(char* MESH,
                       mesh_t &mesh,
                       node_t &node,
                       elem_t &elem,
                       corner_t &corner,
                       const size_t num_dims,
                       const size_t rk_num_bins);

// for string delimiter parsing
std::vector<std::string> split (std::string s, std::string delimiter);

void readVTKPn(char* MESH,
                 mesh_t &mesh,
                 node_t &node,
                 elem_t &elem,
                 zone_t &zone,
                 mat_pt_t &mat_pt,
                 corner_t &corner,
                 fe_ref_elem_t &ref_elem,
                 const size_t num_dims,
               const size_t rk_num_bins);

void input(CArrayKokkos <material_t> &material,
           CArrayKokkos <mat_fill_t> &mat_fill,
           CArrayKokkos <boundary_t> &boundary,
           CArrayKokkos <double> &state_vars,
           size_t &num_materials,
           size_t &num_fills,
           size_t &num_boundaries,
           size_t &num_dims,
           size_t &num_state_vars,
           double &dt_start,
           double &time_final,
           double &dt_max,
           double &dt_min,
           double &dt_cfl,
           double &graphics_dt_ival,
           size_t &graphics_cyc_ival,
           size_t &cycle_stop,
           size_t &rk_num_stages,
           bool &viscosity_cond,
           bool &source_cond);

KOKKOS_FUNCTION
void get_gauss_leg_pt_jacobian(const mesh_t &mesh,
                               const elem_t &elem,
                               const fe_ref_elem_t &ref_elem,
                               const DViewCArrayKokkos <double> &node_coords,
                               DViewCArrayKokkos <double> &gauss_legendre_jacobian,
                               DViewCArrayKokkos <double> &gauss_legendre_det_j,
                               DViewCArrayKokkos <double> &gauss_legendre_jacobian_inverse,
                               const int stage);

KOKKOS_FUNCTION
void get_vol(DViewCArrayKokkos <double> &elem_vol,
             const DViewCArrayKokkos <double> &node_coords,
             const DViewCArrayKokkos <double> &legendre_jacobian_det,
             const mesh_t &mesh,
             const elem_t &elem,
             const fe_ref_elem_t &ref_elem);


void setup(const CArrayKokkos <material_t> &material,
           const CArrayKokkos <mat_fill_t> &mat_fill,
           const CArrayKokkos <boundary_t> &boundary,
           mesh_t &mesh,
           elem_t &elem,
           zone_t &zone,
           mat_pt_t &mat_pt,
           fe_ref_elem_t &ref_elem,
           const DViewCArrayKokkos <double> &node_coords,
           DViewCArrayKokkos <double> &node_vel,
           DViewCArrayKokkos <double> &node_mass,
           const DViewCArrayKokkos <double> &elem_den,
           const DViewCArrayKokkos <double> &elem_pres,
           const DViewCArrayKokkos <double> &elem_stress,
           const DViewCArrayKokkos <double> &elem_sspd,
           const DViewCArrayKokkos <double> &elem_sie,
           const DViewCArrayKokkos <double> &mat_pt_sie,
           const DViewCArrayKokkos <double> &elem_vol,
           const DViewCArrayKokkos <int> &elem_mat_id,
           const DViewCArrayKokkos <double> &elem_statev,
           const CArrayKokkos <double> &state_vars,
           const size_t num_fills,
           const size_t rk_num_bins,
           const size_t num_bdy_sets,
           const size_t num_materials,
           const size_t num_state_vars);

void get_lumped_mass(mesh_t &mesh,
                     fe_ref_elem_t & ref_elem,
                     const DViewCArrayKokkos <double> &DetJac,
                     const DViewCArrayKokkos <double> &den,
                     const DViewCArrayKokkos <double> &M_u,
                     const DViewCArrayKokkos <double> &M_e,
                     DViewCArrayKokkos <double> & nodal_mass,
                     DViewCArrayKokkos <double> & zonal_mass);

void assemble_mass_matrices(const mesh_t &mesh,
                            const fe_ref_elem_t &ref_elem,
                            const DViewCArrayKokkos <double> &rho,
                            const DViewCArrayKokkos <double> &DetJac,
                            DViewCArrayKokkos <double> &M_u,
                            DViewCArrayKokkos <double> &M_e);

void pointwise_mass_conservation(DViewCArrayKokkos <double> &den,
                                 const DViewCArrayKokkos <double> &DetJac,
                                 const DViewCArrayKokkos <double> &den0DetJac0,
                                 const mat_pt_t &mat_pt);

void get_den0DetJac0( DViewCArrayKokkos <double> &den0DetJac0,
                    const DViewCArrayKokkos <double> &den,
                    const DViewCArrayKokkos <double> &DetJac,
                    const mat_pt_t &mat_pt);

void get_energy_residual(const mesh_t &mesh,
                        const DViewCArrayKokkos <double> &M,
                        const DViewCArrayKokkos <double> &zone_sie,
                        const DViewCArrayKokkos <double> &F_e,
                        DViewCArrayKokkos <double> &PSI,
                        const double dt,
                        const int stage,
                        const CArrayKokkos <double> &time_int_weights);

void get_energy_rhs(const mesh_t &mesh,
                    const fe_ref_elem_t &ref_elem,
                    const DViewCArrayKokkos <double> &DetJac,
                    const DViewCArrayKokkos <double> &SigmaJacInv,
                    const DViewCArrayKokkos <double> &node_vel,
                    const DViewCArrayKokkos <double> &stress,
                    const DViewCArrayKokkos <double> &node_coords,
                    DViewCArrayKokkos <double> &F_e,
                    DViewCArrayKokkos <double> &S,
                    const int stage,
                    bool &viscosity_cond,
                    bool &source_cond);

void update_energy(DViewCArrayKokkos <double> &zone_sie,
                   const DViewCArrayKokkos <double> &PSI,
                   const DViewCArrayKokkos <double> &m,
                   const mesh_t &mesh,
                   const int stage);

void get_momentum_rhs(const mesh_t &mesh,
                        const fe_ref_elem_t &ref_elem,
                        const fe_ref_surf_t &ref_surf,
                        const DViewCArrayKokkos <double> &SigmaJacInv,
                        const DViewCArrayKokkos <double> &DetJac,
                        DViewCArrayKokkos <double> &F_u,
                        const int stage,
                        bool &viscosity_cond);

void get_momentum_residual(const mesh_t &mesh,
                            const DViewCArrayKokkos <double> &M,
                            const DViewCArrayKokkos <double> &node_vel,
                            const DViewCArrayKokkos <double> &F_u,
                            DViewCArrayKokkos <double> &PHI,
                            const double dt,
                            const int stage,
                            const CArrayKokkos <double> &time_int_weights);

void update_momentum(DViewCArrayKokkos <double> &node_vel,
                     const DViewCArrayKokkos <double> &PHI,
                     const DViewCArrayKokkos <double> &m,
                     const mesh_t &mesh,
                     const int stage);

void get_stress(const mesh_t &mesh,
                const mat_pt_t &mat_pt,
                const DViewCArrayKokkos <double> &pressure,
                DViewCArrayKokkos <double> &sigma,
                const int stage);

void append_viscosity_to_stress(const mesh_t &mesh,
								const fe_ref_elem_t &ref_elem,
								DViewCArrayKokkos <double> &sigma,
								const DViewCArrayKokkos <double> &node_vel,
                            	const DViewCArrayKokkos <double> &sspd,
                            	const DViewCArrayKokkos <double> &den,
                            	const DViewCArrayKokkos <double> &JacInv,
                            	const DViewCArrayKokkos <double> &Jac,
                            	const DViewCArrayKokkos <double> &DetJac,
                            	const DViewCArrayKokkos <double> &J0Inv,
                            	const DViewCArrayKokkos <double> &h0,
								const int stage);

void get_SigmaJacInv(const mesh_t &mesh,
                     const mat_pt_t &mat_pt,
                     const DViewCArrayKokkos <double> &sigma,
                     const DViewCArrayKokkos <double> &JacInv,
                     DViewCArrayKokkos <double> &SigmaJacInv,
                     const int stage);

void run(const mesh_t &mesh,
         const fe_ref_elem_t &ref_elem,
         const fe_ref_surf_t &ref_surf,
         const elem_t &elem,
         const mat_pt_t &mat_pt,
         CArrayKokkos <boundary_t> &boundary,
         DViewCArrayKokkos <double> &node_vel,
         DViewCArrayKokkos <double> &node_mass,
         DViewCArrayKokkos <double> &M_u,
         DViewCArrayKokkos <double> &zone_sie,
         DViewCArrayKokkos <double> &zone_mass,
         DViewCArrayKokkos <double> &M_e,
         DViewCArrayKokkos <double> &PHI,
         DViewCArrayKokkos <double> &PSI,
         DViewCArrayKokkos <double> &F_u,
         DViewCArrayKokkos <double> &F_e,
         DViewCArrayKokkos <double> &S,
         DViewCArrayKokkos <double> &node_coords,
         DViewCArrayKokkos <double> &Jacobian,
         DViewCArrayKokkos <double> &JacInv,
         DViewCArrayKokkos <double> &DetJac,
         const DViewCArrayKokkos <double> &J0Inv,
         const DViewCArrayKokkos <double> &h0,
         DViewCArrayKokkos <double> &stress,
         DViewCArrayKokkos <double> &SigmaJacInv,
         DViewCArrayKokkos <double> &den,
         DViewCArrayKokkos <double> &den0DetJac0,
         DViewCArrayKokkos <double> &pres,
         DViewCArrayKokkos <double> &sspd,
         DViewCArrayKokkos <double> &elem_vol,
         DViewCArrayKokkos <int> &elem_mat_id,
         const DViewCArrayKokkos <double> &elem_state_vars,
         double &time_value,
         const double time_final,
         const double dt_max,
         const double dt_min,
         const double dt_cfl,
         double &graphics_time,
         size_t graphics_cyc_ival,
         double graphics_dt_ival,
         const size_t cycle_stop,
         const size_t num_stages,
         double dt,
         const double fuzz,
         const double tiny,
         const double small,
         CArray <double> &graphics_times,
         size_t &graphics_id,
         bool &viscosity_cond,
         bool &source_cond);

void write_outputs (const mesh_t &mesh,
                    DViewCArrayKokkos <double> &node_coords,
                    DViewCArrayKokkos <double> &node_vel,
                    DViewCArrayKokkos <double> &node_mass,
                    DViewCArrayKokkos <double> &elem_den,
                    DViewCArrayKokkos <double> &elem_pres,
                    DViewCArrayKokkos <double> &elem_stress,
                    DViewCArrayKokkos <double> &elem_sspd,
                    DViewCArrayKokkos <double> &elem_sie,
                    DViewCArrayKokkos <double> &elem_vol,
                    DViewCArrayKokkos <double> &elem_mass,
                    DViewCArrayKokkos <size_t> &elem_mat_id,
                    CArray <double> &graphics_times, 
                    size_t &graphics_id,
                    const double time_value);


void VTKHexN(const mesh_t &mesh,
             const DViewCArrayKokkos <double> &node_coords,
             const DViewCArrayKokkos <double> &node_vel,
             const DViewCArrayKokkos <double> &node_mass,
             const DViewCArrayKokkos <double> &mat_pt_den,
             const DViewCArrayKokkos <double> &mat_pt_pres,
             const DViewCArrayKokkos <double> &mat_pt_stress,
             const DViewCArrayKokkos <double> &mat_pt_sspd,
             const DViewCArrayKokkos <double> &zone_sie,
             const DViewCArrayKokkos <double> &elem_vol,
             const DViewCArrayKokkos <int> &elem_mat_id,
             CArray <double> &graphics_times,
             size_t &graphics_id,
             const double time_value);


void state_file(const mesh_t &mesh,
                const DViewCArrayKokkos <double> &node_coords,
                const DViewCArrayKokkos <double> &node_vel,
                const DViewCArrayKokkos <double> &mat_pt_h,
                const DViewCArrayKokkos <double> &node_mass,
                const DViewCArrayKokkos <double> &elem_den,
                const DViewCArrayKokkos <double> &elem_pres,
                const DViewCArrayKokkos <double> &elem_stress,
                const DViewCArrayKokkos <double> &elem_sspd,
                const DViewCArrayKokkos <double> &elem_sie,
                const DViewCArrayKokkos <double> &elem_vol,
                const DViewCArrayKokkos <int> &elem_mat_id,
                const double time_value );


void tag_bdys(const CArrayKokkos <boundary_t> &boundary,
              mesh_t &mesh,
              const DViewCArrayKokkos <double> &node_coords);


KOKKOS_FUNCTION
size_t check_bdy(const size_t patch_gid,
                 const int this_bc_tag,
                 const double val,
                 const mesh_t &mesh,
                 const DViewCArrayKokkos <double> &node_coords);


void build_boundry_node_sets(const CArrayKokkos <boundary_t> &boundary,
                             mesh_t &mesh);


void boundary_velocity(const mesh_t &mesh,
                       const CArrayKokkos <boundary_t> &boundary,
                       DViewCArrayKokkos <double> &node_vel,
                       const double time_value);

KOKKOS_FUNCTION 
void eval_sie(const DViewCArrayKokkos <double> sie,
              const int elem_gid,
              const int legendre_lid,
              const mesh_t &mesh,
              const fe_ref_elem_t &ref_elem,
              double &val,
              const int stage);

KOKKOS_FUNCTION 
void eval_vel(const DViewCArrayKokkos <double> vel,
              const int elem_gid,
              const int legendre_lid,
              const mesh_t &mesh,
              const fe_ref_elem_t &ref_elem,
              CArrayKokkos <double> &val,
              const int stage);

void get_J0Inv(const DViewCArrayKokkos <double> &JInv,
               DViewCArrayKokkos <double> &J0Inv,
               const mat_pt_t &mat_pt);

KOKKOS_FUNCTION
void compute_JJ0Inv(const DViewCArrayKokkos <double> &J,
                    const DViewCArrayKokkos <double> &J0Inv,
                    double JJ0Inv[3][3],
                    const int legendre_gid,
                    const mesh_t &mesh);

void add_viscosity_momentum(const mesh_t &mesh,
                            const fe_ref_elem_t &ref_elem,
                            const int stage,
                            const DViewCArrayKokkos <double> &node_vel,
                            const DViewCArrayKokkos <double> &sspd,
                            const DViewCArrayKokkos <double> &den,
                            const DViewCArrayKokkos <double> &JacInv,
                            const DViewCArrayKokkos <double> &Jac,
                            const DViewCArrayKokkos <double> &DetJac,
                            const DViewCArrayKokkos <double> &J0Inv,
                            const DViewCArrayKokkos <double> &h0,
                            DViewCArrayKokkos <double> &F);

void add_viscosity_energy(const mesh_t &mesh,
                            const fe_ref_elem_t &ref_elem,
                            const int stage,
                            const DViewCArrayKokkos <double> &node_vel,
                            const DViewCArrayKokkos <double> &sspd,
                            const DViewCArrayKokkos <double> &den,
                            const DViewCArrayKokkos <double> &JacInv,
                            const DViewCArrayKokkos <double> &Jac,
                            const DViewCArrayKokkos <double> &DetJac,
                            const DViewCArrayKokkos <double> &J0Inv,
                            const DViewCArrayKokkos <double> &h0,
                            DViewCArrayKokkos <double> &F);

KOKKOS_FUNCTION 
void eval_grad_u(const mesh_t &mesh,
                 const fe_ref_elem_t &ref_elem,
                 const int elem_gid,
				 const int leg_lid,
				 const int leg_gid,
                 const DViewCArrayKokkos <double> &node_vel,
                 const DViewCArrayKokkos <double> &JacInv,
                 double grad_u[3][3],
                 const int stage);

KOKKOS_INLINE_FUNCTION
void get_viscosity_coefficient(const mesh_t &mesh,
                               const int leg_gid,
                               double &alpha,
                               const DViewCArrayKokkos <double> &sspd,
                               const DViewCArrayKokkos <double> &den,
                               const DViewCArrayKokkos <double> &J0Inv,
                               const DViewCArrayKokkos <double> &Jac,
                               const double h0,
                               const double grad_u[3][3]);

void get_h0(const DViewCArrayKokkos <double> &elem_vol,
            DViewCArrayKokkos <double> &h0,
            const mesh_t &mesh,
            const fe_ref_elem_t &ref_elem);

void get_vol(DViewCArrayKokkos <double> &elem_vol,
             const DViewCArrayKokkos <double> &node_coords,
             const DViewCArrayKokkos <double> &legendre_jacobian_det,
             const mesh_t &mesh,
             const elem_t &elem,
             const fe_ref_elem_t &ref_elem);

KOKKOS_FUNCTION
void ideal_gas(const DViewCArrayKokkos <double> &elem_pres,
               const DViewCArrayKokkos <double> &elem_stress,
               const size_t elem_gid,
               const size_t legendre_gid,
               const DViewCArrayKokkos <double> &elem_state_vars,
               const DViewCArrayKokkos <double> &elem_sspd,
               const double den,
               const double sie);

void update_thermo(const mesh_t &mesh, 
                   const fe_ref_elem_t &ref_elem,
                   const DViewCArrayKokkos <double> sie,
                   const DViewCArrayKokkos <double> den,
                   DViewCArrayKokkos <double> pres,
                   DViewCArrayKokkos <double> stress,
                   DViewCArrayKokkos <double> sspd,
                   const int stage,
                   const DViewCArrayKokkos <double> &elem_state_vars);

KOKKOS_FUNCTION
void user_eos_model(const DViewCArrayKokkos <double> &elem_pres,
                       const DViewCArrayKokkos <double> &elem_stress,
                       const size_t elem_gid,
                       const size_t legendre_gid,
                       const size_t mat_id,
                       const DViewCArrayKokkos <double> &elem_state_vars,
                       const DViewCArrayKokkos <double> &elem_sspd,
                       const double den,
                       const double sie);


void user_model_init(const DCArrayKokkos <double> &file_state_vars,
                     const size_t num_state_vars,
                     const size_t mat_id,
                     const size_t num_elems);

void update_position_rdh(const int stage,
                         double dt,
                         const mesh_t &mesh,
                         DViewCArrayKokkos <double> &node_coords,
                         const DViewCArrayKokkos <double> &node_vel);



void get_timestep_HexN(const mesh_t &mesh,
                  DViewCArrayKokkos <double> &node_coords,
                  DViewCArrayKokkos <double> &node_vel,
                  DViewCArrayKokkos <double> &mat_pt_sspd,
                  DViewCArrayKokkos <double> &elem_vol,
                  double time_value,
                  const double graphics_time,
                  const double time_final,
                  const double dt_max,
                  const double dt_min,
                  const double dt_cfl,
                  double &dt,
                  const double fuzz);

void init_tn(const mesh_t &mesh,
             DViewCArrayKokkos <double> &node_coords,
             DViewCArrayKokkos <double> &node_vel,
             DViewCArrayKokkos <double> &zone_sie,
             DViewCArrayKokkos <double> &stress,
             const int num_stages);

#endif
