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

#ifndef KINETIC_ENERGY_MINIMIZE_TOPOPT_H
#define KINETIC_ENERGY_MINIMIZE_TOPOPT_H

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
#include "ROL_Elementwise_Reduce.hpp"
#include "FEA_Module_SGH.h"
#include "FEA_Module_Dynamic_Elasticity.h"
#include "Explicit_Solver.h"

class KineticEnergyMinimize_TopOpt : public ROL::Objective<real_t>
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
    FEA_Module_Dynamic_Elasticity* FEM_Dynamic_Elasticity_;
    ROL::Ptr<ROL_MV> ROL_Force;
    ROL::Ptr<ROL_MV> ROL_Velocities;
    ROL::Ptr<ROL_MV> ROL_Gradients;

    bool useLC_; // Use linear form of energy.  Otherwise use quadratic form.

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
    bool   nodal_density_flag_, time_accumulation;
    int    last_comm_step, last_solve_step, current_step;
    size_t nvalid_modules;
    std::vector<FEA_MODULE_TYPE> valid_fea_modules; // modules that may interface with this objective function
    FEA_MODULE_TYPE set_module_type;
    // std::string my_fea_module = "SGH";
    real_t objective_accumulation;

    KineticEnergyMinimize_TopOpt(Explicit_Solver* Explicit_Solver_Pointer, bool nodal_density_flag)
        : useLC_(true)
    {
        Explicit_Solver_Pointer_ = Explicit_Solver_Pointer;

        valid_fea_modules.push_back(FEA_MODULE_TYPE::SGH);
        valid_fea_modules.push_back(FEA_MODULE_TYPE::Dynamic_Elasticity);
        nvalid_modules = valid_fea_modules.size();

        const Simulation_Parameters& simparam = Explicit_Solver_Pointer_->simparam;
        for (const auto& fea_module : Explicit_Solver_Pointer_->fea_modules)
        {
            for (int ivalid = 0; ivalid < nvalid_modules; ivalid++)
            {
                if (fea_module->Module_Type == FEA_MODULE_TYPE::SGH)
                {
                    FEM_SGH_ = dynamic_cast<FEA_Module_SGH*>(fea_module);
                    set_module_type = FEA_MODULE_TYPE::SGH;
                }
                if (fea_module->Module_Type == FEA_MODULE_TYPE::Dynamic_Elasticity)
                {
                    FEM_Dynamic_Elasticity_ = dynamic_cast<FEA_Module_Dynamic_Elasticity*>(fea_module);
                    set_module_type = FEA_MODULE_TYPE::Dynamic_Elasticity;
                }
            }
        }
        nodal_density_flag_ = nodal_density_flag;
        last_comm_step    = last_solve_step = -1;
        current_step      = 0;
        time_accumulation = true;
        objective_accumulation = 0;

        // ROL_Force = ROL::makePtr<ROL_MV>(FEM_->Global_Nodal_Forces);
        if (set_module_type == FEA_MODULE_TYPE::SGH)
        {
            ROL_Velocities = ROL::makePtr<ROL_MV>(FEM_SGH_->node_velocities_distributed);
        }
        if (set_module_type == FEA_MODULE_TYPE::Dynamic_Elasticity)
        {
            ROL_Velocities = ROL::makePtr<ROL_MV>(FEM_Dynamic_Elasticity_->node_velocities_distributed);
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn update
    ///
    /// \brief Evolve the physical simulation
    ///
    /// \param Design vector
    /// \param Update type
    /// \param Iteration
    ///
    /////////////////////////////////////////////////////////////////////////////
    void update(const ROL::Vector<real_t>& z, ROL::UpdateType type, int iter = -1)
    {
        if (set_module_type == FEA_MODULE_TYPE::SGH)
        {
            update_sgh(z, type, iter);
        }
        if (set_module_type == FEA_MODULE_TYPE::Dynamic_Elasticity)
        {
            update_elasticity(z, type, iter);
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn update_elasticity
    ///
    /// \brief Evolve the simulation using dynamic elasticity
    ///
    /// \param Design vector
    /// \param Update type
    /// \param Iteration
    ///
    /////////////////////////////////////////////////////////////////////////////
    void update_elasticity(const ROL::Vector<real_t>& z, ROL::UpdateType type, int iter = -1)
    {
        // debug
        std::ostream& out = std::cout;
        Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

        current_step++;
        ROL::Ptr<const MV>   zp = getVector(z);
        const_host_vec_array design_densities = zp->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);

        if (type == ROL::UpdateType::Initial)
        {
            // This is the first call to update
            // first linear solve was done in FEA class run function already
            FEM_Dynamic_Elasticity_->comm_variables(zp);
            //update deformation variables
            FEM_Dynamic_Elasticity_->update_forward_solve(zp);
            // initial design density data was already communicated for ghost nodes in init_design()
            // decide to output current optimization state
            FEM_Dynamic_Elasticity_->Explicit_Solver_Pointer_->write_outputs();
        }
        else if (type == ROL::UpdateType::Accept)
        {
        }
        else if (type == ROL::UpdateType::Revert)
        {
            // u_ was set to u=S(x) during a trial update
            // and has been rejected as the new iterate
            // Revert to cached value
            // This is a new value of x
            // communicate density variables for ghosts
            FEM_Dynamic_Elasticity_->comm_variables(zp);
            // update deformation variables
            FEM_Dynamic_Elasticity_->update_forward_solve(zp);
            if (Explicit_Solver_Pointer_->myrank == 0)
            {
                *fos << "called Revert" << std::endl;
            }
        }
        else if (type == ROL::UpdateType::Trial)
        {
            // This is a new value of x
            // communicate density variables for ghosts
            FEM_Dynamic_Elasticity_->comm_variables(zp);
            // update deformation variables
            FEM_Dynamic_Elasticity_->update_forward_solve(zp);
            if (Explicit_Solver_Pointer_->myrank == 0)
            {
                *fos << "called Trial" << std::endl;
            }

            // decide to output current optimization state
            FEM_Dynamic_Elasticity_->Explicit_Solver_Pointer_->write_outputs();
        }
        else // ROL::UpdateType::Temp
        // This is a new value of x used for,
        // e.g., finite-difference checks
        {
            if (Explicit_Solver_Pointer_->myrank == 0)
            {
                *fos << "called Temp" << std::endl;
            }
            FEM_Dynamic_Elasticity_->comm_variables(zp);
            FEM_Dynamic_Elasticity_->update_forward_solve(zp);
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn update_sgh
    ///
    /// \brief Evolve the simulation using staggered grid hydrodynamics
    ///
    /// \param Design vector
    /// \param Update type
    /// \param Iteration
    ///
    /////////////////////////////////////////////////////////////////////////////
    void update_sgh(const ROL::Vector<real_t>& z, ROL::UpdateType type, int iter = -1)
    {
        // debug
        std::ostream& out = std::cout;
        Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

        current_step++;
        ROL::Ptr<const MV>   zp = getVector(z);
        const_host_vec_array design_densities = zp->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);

        if (type == ROL::UpdateType::Initial)
        {
            // This is the first call to update
            if(Explicit_Solver_Pointer_->myrank==0)
            {
                *fos << "called SGH Initial" << std::endl;
            }
                
            FEM_SGH_->comm_variables(zp);
            FEM_SGH_->update_forward_solve(zp);
            //initial design density data was already communicated for ghost nodes in init_design()
            //decide to output current optimization state
            //FEM_SGH_->Explicit_Solver_Pointer_->write_outputs();
        }
        else if (type == ROL::UpdateType::Accept)
        {
        }
        else if (type == ROL::UpdateType::Revert)
        {
            // u_ was set to u=S(x) during a trial update
            // and has been rejected as the new iterate
            // Revert to cached value
            // This is a new value of x
            // communicate density variables for ghosts
            if(Explicit_Solver_Pointer_->myrank == 0) *fos << "called SGH Revert" << std::endl;
            
            FEM_SGH_->comm_variables(zp);
            // update deformation variables
            FEM_SGH_->update_forward_solve(zp);
            if (Explicit_Solver_Pointer_->myrank == 0)
            {
                *fos << "called Revert" << std::endl;
            }
        }
        else if (type == ROL::UpdateType::Trial)
        {
            // This is a new value of x
            // communicate density variables for ghosts
            FEM_SGH_->comm_variables(zp);
            // update deformation variables
            FEM_SGH_->update_forward_solve(zp);
            if (Explicit_Solver_Pointer_->myrank == 0)
            {
                *fos << "called Trial" << std::endl;
            }

            // decide to output current optimization state
            // FEM_SGH_->Explicit_Solver_Pointer_->write_outputs();
        }
        else // ROL::UpdateType::Temp
        // This is a new value of x used for,
        // e.g., finite-difference checks
        {
            if (Explicit_Solver_Pointer_->myrank == 0)
            {
                *fos << "called SGH Temp" << std::endl;
            }
            FEM_SGH_->comm_variables(zp);
            FEM_SGH_->update_forward_solve(zp);
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn value
    ///
    /// \brief Returns objective value for optimiation
    ///
    ///
    /// \param Objective value vector
    /// \param Value tolerance
    ///
    /// \return Objective value
    ///
    /////////////////////////////////////////////////////////////////////////////
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

        const_host_vec_array design_densities = zp->getLocalView<HostSpace>(Tpetra::Access::ReadOnly);
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

        // ROL_Force = ROL::makePtr<ROL_MV>(FEM_->Global_Nodal_Forces);
        if (set_module_type == FEA_MODULE_TYPE::SGH)
        {
            ROL_Velocities = ROL::makePtr<ROL_MV>(FEM_SGH_->node_velocities_distributed);
        }
        if (set_module_type == FEA_MODULE_TYPE::Dynamic_Elasticity)
        {
            ROL_Velocities = ROL::makePtr<ROL_MV>(FEM_Dynamic_Elasticity_->node_velocities_distributed);
        }

        std::cout.precision(10);
        if (Explicit_Solver_Pointer_->myrank == 0)
        {
            std::cout << "CURRENT TIME INTEGRAL OF KINETIC ENERGY " << objective_accumulation << std::endl;
        }

        // std::cout << "Ended obj value on task " <<FEM_->myrank  << std::endl;
        return objective_accumulation;
    }

    /////////////////////////////////////////////////////////////////////////////
    ///
    /// \fn gradient
    ///
    /// \brief Calculate design gradient
    ///
    /// \param Vector of gradient values
    /// \param Objective value vector
    /// \param Design tolerance
    ///
    /////////////////////////////////////////////////////////////////////////////
    void gradient(ROL::Vector<real_t>& g, const ROL::Vector<real_t>& z, real_t& tol)
    {
        // std::cout << "Started obj gradient on task " <<FEM_->myrank  << std::endl;
        // get Tpetra multivector pointer from the ROL vector
        ROL::Ptr<const MV> zp = getVector(z);
        ROL::Ptr<MV> gp = getVector(g);

        // communicate ghosts and solve for nodal degrees of freedom as a function of the current design variables
        // FEM_->gradient_print_sync=1;
        // FEM_->gradient_print_sync=0;
        // get local view of the data

        if (set_module_type == FEA_MODULE_TYPE::SGH)
        {
            FEM_SGH_->compute_topology_optimization_gradient_full(zp, gp);
        }
        if (set_module_type == FEA_MODULE_TYPE::Dynamic_Elasticity)
        {
            FEM_Dynamic_Elasticity_->compute_topology_optimization_gradient_full(zp, gp);
        }
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
};

#endif // end header guard
