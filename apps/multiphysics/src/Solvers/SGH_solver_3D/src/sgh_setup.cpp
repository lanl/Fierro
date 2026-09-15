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

#include "sgh_solver_3D.hpp"
//#include "mesh.hpp""
#include "region_fill.hpp"
#include "material.hpp"
#include "boundary_conditions.hpp"
#include "state.hpp"
#include "simulation_parameters.hpp"
#include "geometry_new.hpp"
#include "ELEMENTS.h"




/////////////////////////////////////////////////////////////////////////////
///
/// \fn setup the SGH method
///
/// \brief Allocate state, setup models, and fill mesh regions per the YAML input
///
/////////////////////////////////////////////////////////////////////////////
void SGH3D::setup(SimulationParameters_t& SimulationParamaters, 
                Material_t& Materials, 
                swage::Mesh_t& mesh, 
                BoundaryCondition_t& Boundary,
                State_t& State)
{
    // add a flag on whether SGH was set up, if(SGH_setup_already==false)

    const size_t num_mats = Materials.num_mats; // the number of materials on the mesh

    // Example host-side Logger usage: non-collective, rank-local appends into
    // a buffer that is flushed (collectively) by the Driver at phase
    // boundaries. Safe even inside rank-varying loops because the Logger does
    // NOT touch MPI in the hot path.
    if (log) log->info("Setting up SGH solver, state vars and sspd stress (num_mats=%zu)\n", num_mats);

    // calculate pressure, sound speed, and stress for each material
    for (int mat_id = 0; mat_id < num_mats; mat_id++) {

        // call the initialization function for state vars
        init_state_vars(Materials,
                        mesh,
                        State.MaterialPoints.eos_state_vars,
                        State.MaterialPoints.strength_state_vars,
                        State.MaterialToMeshMaps.elem_in_mat_elem,
                        State.MaterialPoints.num_material_points.host(mat_id),
                        mat_id);

        // call the init function for pressure, sound speed, and stress
        init_press_sspd_stress(Materials,
                               mesh,
                               State.MaterialPoints.den,
                               State.MaterialPoints.pres,
                               State.MaterialPoints.stress,
                               State.MaterialPoints.sspd,
                               State.MaterialPoints.sie,
                               State.MaterialPoints.eos_state_vars,
                               State.MaterialPoints.strength_state_vars,
                               State.MaterialPoints.shear_modulii,
                               State.MaterialPoints.num_material_points.host(mat_id),
                               mat_id);
    } // for loop over mat_id

    // set corner and node masses to zero
    init_corner_node_masses_zero(mesh, State.node.mass, State.corner.mass);

    if (log) log->info("Calculating corner and node masses\n");
    // calculate corner and node masses on the mesh. Each rank sees a
    // rank-local num_mats that may differ; the Logger tolerates this because
    // appends are non-collective.
    for (int mat_id = 0; mat_id < num_mats; mat_id++) {

        if (log) log->info("Calculating corner mass for material %d\n", mat_id);

        // Example of capturing a Logger::Handle by value into a kernel.
        // Obtain `lh` on the host, then pass/capture it into the FOR_ALL.
        // See region_fill.cpp :: log_mat_elem_probe for the actual FOR_ALL
        // using `lh.info(...)` and `FLOG_DEV(lh, INFO, ...)`.
        // if (log && mat_id == 0) {
        //     auto lh = log->handle();
        //     log_mat_elem_probe(lh,
        //                        State.MaterialToMeshMaps.elem_in_mat_elem,
        //                        State.MaterialToMeshMaps.num_mat_elems.host(mat_id),
        //                        static_cast<size_t>(mat_id),
        //                        /* probe_gid = */ 0);
        // }

        calc_corner_mass(Materials,
                         mesh,
                         State.node.coords,
                         State.node.mass,
                         State.corner.mass,
                         State.MaterialPoints.mass,
                         State.MaterialToMeshMaps.elem_in_mat_elem,
                         State.MaterialToMeshMaps.num_mat_elems.host(mat_id),
                         mat_id);
    } // end for mat_id

    if (log) log->info("Calculating node mass\n");

    calc_node_mass(mesh,
                   State.node.coords,
                   State.node.mass,
                   State.corner.mass);

    if (log) log->info("Done calculating node mass\n");

    // std::cout << "Setting up fracture" << std::endl;
    // setting up fracture
    for (size_t i = 0; i < mesh.num_bdy_sets; i++) {
        // if fracture is allowed, then set up the fracture bank
        // note, allow_fracture is set in the parse_bdy_conds_inputs.cpp file and boundary_conditions.h file
        // checking if fracture is allowed... if = 0 then fracture is not enabled; if = 1, then fracture is enabled:
        // printf("Boundary.allow_fracture = %d\n", Boundary.allow_fracture);
        if (Boundary.allow_fracture) {
            // printf("Setting up global fracture (cohesive zones)\n");
            doing_fracture = true;
        
            // calling initialize for the cohesive zones bank
            //printf("Calling initialize()...\n");
            //cohesive_zones_t cohesive_zones_bank;

            // this->cohesive_zones_bank.initialize(mesh, State, SimulationParamaters);
            this->cohesive_zones_bank.initialize(
                mesh.nodes_in_elem,
                mesh.elems_in_node,
                mesh.bdy_nodes,
                mesh.num_bdy_nodes,
                State.node.coords,
                cohesive_zones_bank.geom_tol
            );

            // done calling initialize
            // printf("Done calling initialize()...\n");
            break; 
        }
    }
    // end setting up fracture
    
    // Setting up contact
    if (log) log->info("Setting up contact\n");
    // todo: should this be handled inside of src/boundary_conditions/stress/global_contact ?
    for (size_t i = 0; i < mesh.num_bdy_sets; i++) {
        if (Boundary.allow_preload) {
            if (log) log->info("Setting up preload contact\n");
            doing_preload = true;
            doing_contact = true;
            break;
        }
        if (Boundary.allow_contact) {
            if (log) log->info("Setting up global contact\n");
            doing_contact = true;
            break;
        }
    }


    // Setup the reference element
       // ================================================================
    // Create quadrature along with the reference element and surface

    if (log) log->info("Building reference elements and quadrature in SGH setup\n");

    // the minimum quadrature for FE hydrodynamics based on elem order
    const size_t num_DOFs_1d = 2; // Limited to linear elements for SGH solver
    const size_t num_qpts_1d = 2; // hard coded for SGH solver
    const size_t elem_dims = 3; // 3D elements
    const size_t elem_order = 1; // linear elements


    // ---- reference element ----

    // create quadrature
    Quad.initialize_quadrature(reference_space::GaussLegendre,
                               num_qpts_1d,
                               elem_dims);

    // p_order is the basis order for the Lagrange polynomial defining the element
    FERefElem.initialize_ref_elem(reference_space::arbitraryOrderElement,
                                  reference_space::LagrangeLobatto,
                                  Quad,
                                  elem_order);    

    // ---- reference surface ----
    SurfQuad.initialize_quadrature(reference_space::GaussLegendre, 
                                   num_qpts_1d, 
                                   elem_dims); 

    RefSurf.initialize_ref_surf(SurfQuad,
                                FERefElem);

    // Map to get from quadrature points on the surface to the element
    int num_surfaces = mesh.num_surfs;
    const size_t num_surf_qpts = SurfQuad.num_qpts_in_surf;
    this->surf_qpt_qpt_map = CArrayKokkos<int>(num_surfaces, 2, num_surf_qpts, "surf_qpt_qpt_map");
    this->surf_qpt_qpt_map.set_values(-1);

    build_quadrature_point_connectivity(mesh, RefSurf, this->surf_qpt_qpt_map, State.node.coords); 



    // Setup basis tables
    if (log) log->info("Setting up basis tables\n");
    const size_t num_nodes_in_elem = mesh.num_nodes_in_elem;
    const size_t num_qpts_in_elem = Quad.num_qpts_in_elem;
    const size_t num_surfs_in_elem = mesh.num_surfs_in_elem;
    const size_t num_qpts_in_surf = SurfQuad.num_qpts_in_surf;

    // Allocate the ALE specific jacobian, inverse jacobian, and RHS fluxes and fields
    this->elem_det_jac = CArrayKokkos<double>(mesh.num_elems, num_qpts_in_elem, "elem_det_jacobian");
    this->inv_jac_ijq = CArrayKokkos<double>(mesh.num_elems, elem_dims, elem_dims, num_qpts_in_elem, "inv_jac_ijq");
    
    this->surf_vn = CArrayKokkos<double>(mesh.num_surfs, num_qpts_in_surf, "surf_vn");
    this->RHS_surf_flux = DRaggedRightArrayKokkos<double>(State.MaterialToMeshMaps.num_mat_elems_buffer, num_surfs_in_elem, num_qpts_in_surf, "RHS_surf_flux");
    this->RHS_corner = DRaggedRightArrayKokkos<double>(State.MaterialCorners.num_material_corners_buffer, "RHS_corner");

    this->qpt_adv_vel = CArrayKokkos<double>(mesh.num_elems, num_qpts_in_elem, elem_dims, "qpt_adv_vel");
    this->mat_qpt_field = DRaggedRightArrayKokkos<double>(State.MaterialToMeshMaps.num_mat_elems_buffer, num_qpts_in_elem, "mat_qpt_field");

    // Setup the basis tables
    tables.basis_row_sum = CArrayKokkos<double>(num_qpts_in_elem, "basis_row_sum");
    
    
    FOR_ALL(qpt_lid, 0, num_qpts_in_elem, {
        double sum = 0.0;
        for(size_t dof_lid = 0; dof_lid < num_nodes_in_elem; dof_lid++){
            sum += FERefElem.qpt_basis(qpt_lid, dof_lid);
        }
        tables.basis_row_sum(qpt_lid) = sum;
    });
    
    tables.basis_dq       = CArrayKokkos<double>(num_nodes_in_elem, num_qpts_in_elem, "basis_dq");
    tables.grad_basis_njq = CArrayKokkos<double>(num_nodes_in_elem, elem_dims, num_qpts_in_elem, "grad_basis_njq");
    tables.grad_basis_qjd = CArrayKokkos<double>(num_qpts_in_elem, elem_dims, num_nodes_in_elem, "grad_basis_qjd");
    
    FOR_ALL(qpt_lid, 0, num_qpts_in_elem, {
        for(size_t dof_lid = 0; dof_lid < num_nodes_in_elem; dof_lid++){
            tables.basis_dq(dof_lid, qpt_lid) = FERefElem.qpt_basis(qpt_lid, dof_lid);
            for(size_t dim = 0; dim < elem_dims; dim++){
                const double val = FERefElem.qpt_grad_basis(qpt_lid, dof_lid, dim);
                tables.grad_basis_njq(dof_lid, dim, qpt_lid) = val;
                tables.grad_basis_qjd(qpt_lid, dim, dof_lid) = val;
            }
        }
    });
    tables.surf_basis_fdq = CArrayKokkos<double>(num_surfs_in_elem, num_nodes_in_elem,
                                                 num_qpts_in_surf, "surf_basis_fdq");
    tables.surf_grad_fjdq = CArrayKokkos<double>(num_surfs_in_elem, elem_dims, num_nodes_in_elem,
                                                 num_qpts_in_surf, "surf_grad_fjdq");
    FOR_ALL(face_lid, 0, num_surfs_in_elem, {
        for(size_t qpt_lid = 0; qpt_lid < num_qpts_in_surf; qpt_lid++){
            for(size_t dof_lid = 0; dof_lid < num_nodes_in_elem; dof_lid++){
                tables.surf_basis_fdq(face_lid, dof_lid, qpt_lid) =
                    RefSurf.qpt_basis(face_lid, qpt_lid, dof_lid);
                for(size_t dim = 0; dim < elem_dims; dim++){
                    tables.surf_grad_fjdq(face_lid, dim, dof_lid, qpt_lid) =
                        RefSurf.qpt_grad_basis(face_lid, qpt_lid, dof_lid, dim);
                }
            }
        }
    });
    Kokkos::fence();

    // ================================================================
    // Step 1: build the volume matrix for nodal DG at t=0

    build_element_geometry(mesh, tables, State.node.coords, elem_det_jac, inv_jac_ijq,
        mesh.num_elems, num_qpts_in_elem, num_nodes_in_elem);
    build_lumped_volume(mesh, FERefElem, Quad, tables, elem_det_jac, State.corner.volume,
        mesh.num_elems, num_qpts_in_elem, num_nodes_in_elem);
    Kokkos::fence();

    this->mesh_node_target_coords = CArrayKokkos<double>(mesh.num_nodes, mesh.num_dims, "mesh_node_target_coords");

    FOR_ALL(node_gid, 0, mesh.num_nodes, {
        for(size_t dim = 0; dim < mesh.num_dims; dim++){
            this->mesh_node_target_coords(node_gid, dim) = State.node.coords(node_gid, dim);
        }
    });

    // Initialize the mesh node velocity with the same communication plan as the node velocity
    this->mesh_node_velocity = MPICArrayKokkos<double>(mesh.num_nodes, mesh.num_dims, "mesh_node_velocity");
    if (State.node.vel.comm_plan_ != nullptr) {
        this->mesh_node_velocity.initialize_comm_plan(*State.node.vel.comm_plan_);
    }

    
} // end SGH setup
