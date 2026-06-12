#include "ao_sgh_solver_3D.hpp"
#include "boundary_conditions.hpp"

namespace ao_sgh
{

void AO_SGH3D::boundary_velocity(const swage::Mesh&         mesh,
                                 const BoundaryCondition_t& BoundaryConditions,
                                 DCArrayKokkos<double>&     node_vel,
                                 const double               time_value) const
{
    const size_t num_vel_bdy_sets =
        BoundaryConditions.num_vel_bdy_sets_in_solver.host(this->solver_id);

    for (size_t bc_lid = 0; bc_lid < num_vel_bdy_sets; bc_lid++) {

        const size_t bdy_set =
            BoundaryConditions.vel_bdy_sets_in_solver.host(this->solver_id, bc_lid);

        FOR_ALL(bdy_node_lid, 0, mesh.num_bdy_nodes_in_set.host(bdy_set), {
            const size_t bdy_node_gid =
                mesh.bdy_nodes_in_set(bdy_set, bdy_node_lid);

            BoundaryConditions.BoundaryConditionFunctions(bdy_set).velocity(
                mesh,
                BoundaryConditions.BoundaryConditionEnums,
                BoundaryConditions.velocity_bc_global_vars,
                BoundaryConditions.bc_state_vars,
                node_vel,
                time_value,
                1,
                bdy_node_gid,
                bdy_set);
        });
    }
}

void AO_SGH3D::boundary_force(const swage::Mesh&         mesh,
                              const BoundaryCondition_t& BoundaryConditions,
                              DCArrayKokkos<double>&     node_force) const
{
    const size_t num_vel_bdy_sets =
        BoundaryConditions.num_vel_bdy_sets_in_solver.host(this->solver_id);

    for (size_t bc_lid = 0; bc_lid < num_vel_bdy_sets; bc_lid++) {

        const size_t bdy_set =
            BoundaryConditions.vel_bdy_sets_in_solver.host(this->solver_id, bc_lid);

        FOR_ALL(bdy_node_lid, 0, mesh.num_bdy_nodes_in_set.host(bdy_set), {
            const size_t bdy_node_gid =
                mesh.bdy_nodes_in_set(bdy_set, bdy_node_lid);

            double mag2 = 0.0;
            for (size_t dim = 0; dim < mesh.num_dims; dim++) {
                const double nd = BoundaryConditions.velocity_bc_global_vars(bdy_set, dim);
                mag2 += nd * nd;
            }
            if (mag2 == 0.0) return;
            // Project out the component along the BC direction:
            //   f := f - (f.n / |n|^2) * n.
            double fdotn = 0.0;
            for (size_t dim = 0; dim < mesh.num_dims; dim++) {
                const double nd = BoundaryConditions.velocity_bc_global_vars(bdy_set, dim);
                fdotn += node_force(bdy_node_gid, dim) * nd;
            }
            const double s = fdotn / mag2;
            for (size_t dim = 0; dim < mesh.num_dims; dim++) {
                const double nd = BoundaryConditions.velocity_bc_global_vars(bdy_set, dim);
                node_force(bdy_node_gid, dim) -= s * nd;
            }
        });
    }
}

} // namespace ao_sgh
