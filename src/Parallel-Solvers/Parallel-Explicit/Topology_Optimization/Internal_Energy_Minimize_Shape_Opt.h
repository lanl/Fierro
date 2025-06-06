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

#ifndef INTERNAL_ENERGY_MINIMIZE_SHAPE_OPT_H
#define INTERNAL_ENERGY_MINIMIZE_SHAPE_OPT_H

#include "matar.h"
#include "elements.h"
#include <string>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"

#include "ROL_Types.hpp"
#include <ROL_TpetraMultiVector.hpp>
#include "ROL_Objective.hpp"
#include "Fierro_Optimization_Objective.hpp"
#include "ROL_Elementwise_Reduce.hpp"
#include "FEA_Module_SGH.h"
#include "Explicit_Solver.h"

class InternalEnergyMinimize_ShapeOpt : public FierroOptimizationObjective
{
typedef Tpetra::Map<>::local_ordinal_type LO;
typedef Tpetra::Map<>::global_ordinal_type GO;
typedef Tpetra::Map<>::node_type Node;
typedef Tpetra::Map<LO, GO, Node> Map;
typedef Tpetra::MultiVector<real_t, LO, GO, Node> MV;
typedef ROL::Vector<real_t> V;
typedef ROL::TpetraMultiVector<real_t, LO, GO, Node> ROL_MV;

using traits = Kokkos::ViewTraits<LO*, Kokkos::LayoutLeft, void, void>;
using array_layout    = typename traits::array_layout;
using execution_space = typename traits::execution_space;
using device_type     = typename traits::device_type;
using memory_traits   = typename traits::memory_traits;
using global_size_t   = Tpetra::global_size_t;

typedef Kokkos::View<real_t*, Kokkos::LayoutRight, device_type, memory_traits> values_array;
typedef Kokkos::View<GO*, array_layout, device_type, memory_traits> global_indices_array;
typedef Kokkos::View<LO*, array_layout, device_type, memory_traits> indices_array;

// typedef Kokkos::DualView<real_t**, Kokkos::LayoutLeft, device_type>::t_dev vec_array;
typedef MV::dual_view_type::t_dev vec_array;
typedef MV::dual_view_type::t_host host_vec_array;
typedef Kokkos::View<const real_t**, array_layout, HostSpace, memory_traits> const_host_vec_array;
typedef Kokkos::View<const real_t**, array_layout, device_type, memory_traits> const_vec_array;
typedef MV::dual_view_type dual_vec_array;

private:

    Explicit_Solver* Explicit_Solver_Pointer_;
    FEA_Module_SGH*  FEM_SGH_;
    ROL::Ptr<ROL_MV> ROL_Force;
    ROL::Ptr<ROL_MV> ROL_Velocities;
    ROL::Ptr<ROL_MV> ROL_Gradients;
    Teuchos::RCP<MV> previous_gradients;
    real_t initial_internal_energy;
    real_t previous_objective_accumulation, objective_sign;

    bool first_init; //prevents ROL from calling init computation twice at start for the AL algorithm

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn getVector
    ///
    /// \brief Retrieves ROL vector at desired location
    ///
    /// \param Pointer to desired ROL vector
    ///
    /// \return Returns ROL MV vector
    ///
    /////////////////////////////////////////////////////////////////////////////
    ROL::Ptr<const MV> getVector(const V& x)
    {
        return dynamic_cast<const ROL_MV&>(x).getVector();
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn getVector
    ///
    /// \brief Retrieves ROL vector at desired location
    ///
    /// \param Pointer to desired ROL vector
    ///
    /// \return Returns const ROL MV vector
    ///
    /////////////////////////////////////////////////////////////////////////////
    ROL::Ptr<MV> getVector(V& x)
    {
        return dynamic_cast<ROL_MV&>(x).getVector();
    }

public:
    int last_comm_step, last_solve_step, current_step;
    int num_dim;
    size_t nvalid_modules;
    size_t nlocal_nodes, num_corners, num_nodes_in_elem, rnum_elem;
    DViewCArrayKokkos<double> node_mass, node_coords;
    std::vector<FEA_MODULE_TYPE> valid_fea_modules; // modules that may interface with this objective function
    FEA_MODULE_TYPE set_module_type;
    // std::string my_fea_module = "SGH";

    InternalEnergyMinimize_ShapeOpt(Explicit_Solver* Explicit_Solver_Pointer)
    {
        Explicit_Solver_Pointer_ = Explicit_Solver_Pointer;
        first_init = false;
        valid_fea_modules.push_back(FEA_MODULE_TYPE::SGH);
        nvalid_modules = valid_fea_modules.size();
        objective_sign = 1;
        num_dim = Explicit_Solver_Pointer_->simparam.num_dims;
        const Simulation_Parameters& simparam = Explicit_Solver_Pointer_->simparam;
        for (const auto& fea_module : Explicit_Solver_Pointer_->fea_modules) {
            for (int ivalid = 0; ivalid < nvalid_modules; ivalid++) {
                if (fea_module->Module_Type == FEA_MODULE_TYPE::SGH) {
                    FEM_SGH_ = dynamic_cast<FEA_Module_SGH*>(fea_module);
                    set_module_type = FEA_MODULE_TYPE::SGH;
                    node_mass = FEM_SGH_->node_mass;
                    node_coords = FEM_SGH_->node_coords;
                    nlocal_nodes = FEM_SGH_->nlocal_nodes;
                    num_corners = FEM_SGH_->num_corners;
                    num_nodes_in_elem = FEM_SGH_->num_nodes_in_elem;
                    rnum_elem = FEM_SGH_->rnum_elem;
                }
            }
        }
        last_comm_step    = last_solve_step = -1;
        current_step      = 0;
        time_accumulation = true;
        
        previous_gradients = Teuchos::rcp(new MV(Explicit_Solver_Pointer_->map, num_dim));
        if(Explicit_Solver_Pointer_->simparam.optimization_options.maximize_flag){
            objective_sign = -1;
        }
        // ROL_Force = ROL::makePtr<ROL_MV>(FEM_->Global_Nodal_Forces);
        if (set_module_type == FEA_MODULE_TYPE::SGH) {
            ROL_Velocities = ROL::makePtr<ROL_MV>(FEM_SGH_->node_velocities_distributed);
        }
    }

  /* --------------------------------------------------------------------------------------
   Update solver state variables to synchronize with the current design variable vector, z
  ----------------------------------------------------------------------------------------- */

    void update(const ROL::Vector<real_t>& z, ROL::UpdateType type, int iter = -1)
    {
        if (set_module_type == FEA_MODULE_TYPE::SGH) {
            update_sgh(z, type, iter);
        }
    }

  /* --------------------------------------------------------------------------------------
   Update solver state variables to synchronize with the current design variable vector, z
  ----------------------------------------------------------------------------------------- */

    void update_sgh(const ROL::Vector<real_t>& z, ROL::UpdateType type, int iter = -1)
    {
        // debug
        std::ostream& out = std::cout;
        Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
        bool print_flag = false;
        ROL::Ptr<const MV>   zp = getVector(z); //tpetra multivector wrapper on design vector
        const_host_vec_array design_coordinates = zp->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);

        if (type == ROL::UpdateType::Initial) {
            if(first_init){
                // This is the first call to update
                if (Explicit_Solver_Pointer_->myrank == 0) {
                    *fos << "called SGH Initial" << std::endl;
                }

                FEM_SGH_->comm_variables(zp);
                FEM_SGH_->update_forward_solve(zp);
                if(Explicit_Solver_Pointer_->myrank == 0){
                    std::cout << "CURRENT TIME INTEGRAL OF INTERNAL ENERGY " << objective_accumulation << std::endl;
                }
                FEM_SGH_->compute_shape_optimization_adjoint_full(zp);
                previous_objective_accumulation = objective_accumulation;
                previous_gradients->assign(*(FEM_SGH_->cached_design_gradients_distributed));
                // initial design density data was already communicated for ghost nodes in init_design()
                // decide to output current optimization state
                // FEM_SGH_->Explicit_Solver_Pointer_->write_outputs();
            }
            first_init = true;
        }
        else if (type == ROL::UpdateType::Accept) {
            if (Explicit_Solver_Pointer_->myrank == 0) {
                *fos << "called Accept" << std::endl;
            }
            
            previous_objective_accumulation = objective_accumulation;
            previous_gradients->assign(*(FEM_SGH_->cached_design_gradients_distributed));
        }
        else if (type == ROL::UpdateType::Revert) {
            // u_ was set to u=S(x) during a trial update
            // and has been rejected as the new iterate
            // Revert to cached value
            // This is a new value of x
            // communicate density variables for ghosts
            if (Explicit_Solver_Pointer_->myrank == 0) { *fos << "called Revert" << std::endl; }
            objective_accumulation = previous_objective_accumulation;
            FEM_SGH_->cached_design_gradients_distributed->assign(*previous_gradients);
            if(Explicit_Solver_Pointer_->myrank == 0){
                std::cout << "CURRENT TIME INTEGRAL OF INTERNAL ENERGY " << objective_accumulation << std::endl;
            }
            // FEM_SGH_->comm_variables(zp);
            // // update deformation variables
            // FEM_SGH_->update_forward_solve(zp);
            // FEM_SGH_->compute_shape_optimization_adjoint_full(zp);
        }
        else if (type == ROL::UpdateType::Trial) {
            // This is a new value of x
            current_step++;
            if(current_step%FEM_SGH_->simparam->optimization_options.optimization_output_freq==0){
                print_flag = true;
            }
            if (Explicit_Solver_Pointer_->myrank == 0) {
                *fos << "called Trial" << std::endl;
            }
            // communicate density variables for ghosts
            FEM_SGH_->comm_variables(zp);
            // update deformation variables
            FEM_SGH_->update_forward_solve(zp, print_flag);
            
            if(Explicit_Solver_Pointer_->myrank == 0){
                std::cout << "CURRENT TIME INTEGRAL OF INTERNAL ENERGY " << objective_accumulation << std::endl;
            }
            FEM_SGH_->compute_shape_optimization_adjoint_full(zp);
            // decide to output current optimization state
            // FEM_SGH_->Explicit_Solver_Pointer_->write_outputs();
        }
        else{ // ROL::UpdateType::Temp
            // This is a new value of x used for,
            // e.g., finite-difference checks
            if (Explicit_Solver_Pointer_->myrank == 0) {
                *fos << "called Temp" << std::endl;
            }
            FEM_SGH_->comm_variables(zp);
            FEM_SGH_->update_forward_solve(zp);
            if(Explicit_Solver_Pointer_->myrank == 0){
                std::cout << "CURRENT TIME INTEGRAL OF INTERNAL ENERGY " << objective_accumulation << std::endl;
            }
            FEM_SGH_->compute_shape_optimization_adjoint_full(zp);
        }
    }

  /* --------------------------------------------------------------------------------------
   Compute time integral contribution for this objective function form
  ----------------------------------------------------------------------------------------- */
    void step_accumulation(const real_t& dt, const size_t& cycle, const size_t& rk_level) {
        
        const_vec_array current_elem_sie;
        const_vec_array previous_elem_sie;
        size_t nlocal_elem_non_overlapping = FEM_SGH_->Explicit_Solver_Pointer_->nlocal_elem_non_overlapping;
        auto nodes_in_elem = FEM_SGH_->nodes_in_elem;
        auto elem_mass = FEM_SGH_->elem_mass;
        bool use_solve_checkpoints    = FEM_SGH_->simparam->optimization_options.use_solve_checkpoints;
        if(use_solve_checkpoints){
            current_elem_sie = FEM_SGH_->element_internal_energy_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            previous_elem_sie = FEM_SGH_->previous_element_internal_energy_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
        }
        else{
            auto forward_solve_internal_energy_data = FEM_SGH_->forward_solve_internal_energy_data;
            current_elem_sie = (*forward_solve_internal_energy_data)[cycle + 1]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            previous_elem_sie = (*forward_solve_internal_energy_data)[cycle]->getLocalView<device_type>(Tpetra::Access::ReadOnly);
        }

        double IE_sum = 0.0;
        double IE_loc_sum = 0.0;
        // extensive KE
        if(FEM_SGH_->simparam->optimization_options.optimization_objective_regions.size()){
            int nobj_volumes = FEM_SGH_->simparam->optimization_options.optimization_objective_regions.size();
            auto optimization_objective_regions = FEM_SGH_->simparam->optimization_options.optimization_objective_regions;
            const_vec_array all_initial_node_coords = FEM_SGH_->all_initial_node_coords_distributed->getLocalView<device_type>(Tpetra::Access::ReadOnly);
            FOR_REDUCE_SUM_CLASS(elem_gid, 0, nlocal_elem_non_overlapping, IE_loc_sum, {
                int node_id;
                double current_elem_coords[3];
                bool contained = false;
                current_elem_coords[0] = 0;
                current_elem_coords[1] = 0;
                current_elem_coords[2] = 0;
                for(int inode=0; inode< num_nodes_in_elem; inode++){
                    node_id = nodes_in_elem(elem_gid, inode);
                    current_elem_coords[0] += all_initial_node_coords(node_id, 0)/num_nodes_in_elem;
                    current_elem_coords[1] += all_initial_node_coords(node_id, 1)/num_nodes_in_elem;
                    current_elem_coords[2] += all_initial_node_coords(node_id, 2)/num_nodes_in_elem;
                }
                for(int ivolume = 0; ivolume < nobj_volumes; ivolume++){
                    if(optimization_objective_regions(ivolume).contains(current_elem_coords)){
                        contained = true;
                    }
                }
                if(contained){
                    IE_loc_sum += elem_mass(elem_gid) *  0.5 * (current_elem_sie(elem_gid,0)+previous_elem_sie(elem_gid,0));
                }

                
            }, IE_sum);
        }
        else{
            FOR_REDUCE_SUM_CLASS(elem_gid, 0, nlocal_elem_non_overlapping, IE_loc_sum, {
                IE_loc_sum += elem_mass(elem_gid) * 0.5 * (current_elem_sie(elem_gid,0)+previous_elem_sie(elem_gid,0));
            }, IE_sum);
        }

        Kokkos::fence();
        objective_accumulation += IE_sum * dt;
    }

  /* --------------------------------------------------------------------------------------
   Update objective value with the current design variable vector, z
  ----------------------------------------------------------------------------------------- */

    real_t value(const ROL::Vector<real_t>& z, real_t& tol)
    {
        // std::cout << "Started obj value on task " <<FEM_->myrank  << std::endl;
        ROL::Ptr<const MV> zp = getVector(z);
        real_t c = 0.0;

        // debug print
        // std::ostream &out = std::cout;
        // Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
        // if(FEM_->myrank==0)
        // *fos << "Value function z:" << std::endl;
        // zp->describe(*fos,Teuchos::VERB_EXTREME);
        // *fos << std::endl;
        // std::fflush(stdout);

        const_host_vec_array design_coordinates = zp->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
        // communicate ghosts and solve for nodal degrees of freedom as a function of the current design variables
        /*
        if(last_comm_step!=current_step){
          FEM_->comm_variables(zp);
          last_comm_step = current_step;
        }

        if(last_solve_step!=current_step){
          //std::cout << "UPDATED velocities" << std::endl;
          FEM_->update_linear_solve(zp);
          last_solve_step = current_step;
        }
        */
        // debug print of velocities
        // std::ostream &out = std::cout;
        // Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
        // if(FEM_->myrank==0)
        // *fos << "Displacement data :" << std::endl;
        // FEM_->node_velocities_distributed->describe(*fos,Teuchos::VERB_EXTREME);
        // *fos << std::endl;
        // std::fflush(stdout);

        std::cout.precision(10);

        // std::cout << "Ended obj value on task " <<FEM_->myrank  << std::endl;
        return objective_sign*objective_accumulation;
    }

  /* --------------------------------------------------------------------------------------
   Update gradient vector (g) with the current design variable vector, z
  ----------------------------------------------------------------------------------------- */

    void gradient(ROL::Vector<real_t>& g, const ROL::Vector<real_t>& z, real_t& tol)
    {
        // std::cout << "Started obj gradient on task " <<FEM_->myrank  << std::endl;
        // get Tpetra multivector pointer from the ROL vector
        ROL::Ptr<const MV> zp = getVector(z); //pointer to design vector
        ROL::Ptr<MV> gp = getVector(g); //pointer to gradient vector

        // communicate ghosts and solve for nodal degrees of freedom as a function of the current design variables
        // FEM_->gradient_print_sync=1;
        // FEM_->gradient_print_sync=0;
        // get local view of the data

        if (set_module_type == FEA_MODULE_TYPE::SGH) {
            FEM_SGH_->compute_shape_optimization_gradient_full(zp, gp);
        }
        gp->scale(objective_sign);
        // debug print of gradient
        // std::ostream &out = std::cout;
        // Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
        // if(FEM_->myrank==0)
        // *fos << "Gradient data :" << std::endl;
        // gp->describe(*fos,Teuchos::VERB_EXTREME);
        // *fos << std::endl;
        // std::fflush(stdout);
        // for(int i = 0; i < FEM_->nlocal_nodes; i++){
        // objective_gradients(i,0) *= -1;
        // }

        // std::cout << "Objective Gradient called"<< std::endl;
        // debug print of design variables
        // std::ostream &out = std::cout;
        // Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
        // if(FEM_->myrank==0)
        // *fos << "Gradient data :" << std::endl;
        // gp->describe(*fos,Teuchos::VERB_EXTREME);

        // *fos << std::endl;
        // std::fflush(stdout);
        // std::cout << "ended obj gradient on task " <<FEM_->myrank  << std::endl;
    }
    
  //contributes to rate of change of adjoint vector due to term with velocity gradient of objective
    void sie_gradient_adjoint_contribution(vec_array& adjoint_rate_vector, const DViewCArrayKokkos<double>& node_mass,
                                                const DViewCArrayKokkos<double>& elem_mass, const DViewCArrayKokkos<double>& node_vel,
                                                const DViewCArrayKokkos<double>& node_coords, const DViewCArrayKokkos<double>& elem_sie,
                                                const size_t& rk_level){
        

        FOR_ALL_CLASS(elem_gid, 0, rnum_elem, {
            adjoint_rate_vector(elem_gid, 0) = elem_mass(elem_gid);
        }); // end parallel for
        Kokkos::fence();
    }

  //contribution to gradient tally from objective w.r.t design coordinate
    void design_coordinate_gradient_term(vec_array& gradient_vector, const DViewCArrayKokkos<double>& node_mass,
                               const DViewCArrayKokkos<double>& elem_mass, const DViewCArrayKokkos<double>& node_vel,
                               const DViewCArrayKokkos<double>& node_coords, const DViewCArrayKokkos<double>& elem_sie,
                               const size_t& rk_level, const real_t& global_dt = 0){
        auto optimization_objective_regions = FEM_SGH_->simparam->optimization_options.optimization_objective_regions;
        auto nodes_in_elem = FEM_SGH_->nodes_in_elem;
        auto corner_value_storage = FEM_SGH_->corner_value_storage;
        auto corners_in_node = FEM_SGH_->corners_in_node;
        auto num_corners_in_node = FEM_SGH_->num_corners_in_node;
    }

};

#endif // end header guard
