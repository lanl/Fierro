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
 
#ifndef IMPLICIT_SOLVER_H
#define IMPLICIT_SOLVER_H

#include "Solver.h"
#include "Tpetra_computeRowAndColumnOneNorms_decl.hpp"
#include "Tpetra_Details_EquilibrationInfo.hpp"

//#include <Xpetra_Operator.hpp>
//#include <MueLu.hpp>

//forward declarations
namespace MueLu{
  template<class floattype, class local_ind, class global_ind, class nodetype> 
  class Hierarchy;
}

namespace Xpetra{
  template<class floattype, class local_ind, class global_ind, class nodetype> 
  class Operator;
  
  template<class floattype, class local_ind, class global_ind, class nodetype> 
  class MultiVector;

  template<class floattype, class local_ind, class global_ind, class nodetype> 
  class Matrix;
}

class FEA_Module;

class Implicit_Solver: public Solver{

public:
  Implicit_Solver();
  ~Implicit_Solver();

  void run(int argc, char *argv[]);

  //void read_mesh_ensight(char *MESH, bool convert_node_order);

  //void read_mesh_tecplot(char *MESH);

  void read_mesh_ansys_dat(char *MESH);
  
  //setup ghosts and element maps
  void init_maps();

  //void repartition_nodes();

  void init_design();

  void collect_information();

  void sort_information();

  //process input to decide TO problem and FEA modules
  void FEA_module_setup();

  void setup_optimization_problem();
  
  //initialize data for boundaries of the model and storage for boundary conditions and applied loads
  void init_boundaries();
  
  //interfaces between user input and creating data structures for bcs
  void topology_conditions();

  //void vtk_writer();

  //void ensight_writer();

  //interfaces between user input and creating data structures for topology conditions
  void generate_tcs();

  void init_topology_conditions (int num_sets);

  void tecplot_writer();

  void parallel_tecplot_writer();

  //void init_boundary_sets(int num_boundary_sets);

  void tag_boundaries(int this_bc_tag, real_t val, int bdy_set, real_t *patch_limits = NULL);

  int check_boundary(Node_Combination &Patch_Nodes, int this_bc_tag, real_t val, real_t *patch_limits);

  //output stream
  Teuchos::RCP<Teuchos::FancyOStream> fos;
  
  swage::mesh_t *init_mesh;
  swage::mesh_t *mesh;
  

  //class Simulation_Parameters *simparam;
  class Simulation_Parameters_Topology_Optimization *simparam_TO;

  //set of enabled FEA modules
  std::vector<std::string> fea_module_types;
  std::vector<FEA_Module*> fea_modules;
  std::vector<bool> fea_module_must_read;
  int nfea_modules;
  int displacement_module;
  
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> Topology_Condition_Patches; //set of patches corresponding to each boundary condition
  CArrayKokkos<size_t, array_layout, device_type, memory_traits> NTopology_Condition_Patches;

  //types of boundary conditions
  enum tc_type {NONE, TO_SURFACE_CONSTRAINT, TO_BODY_CONSTRAINT};

  //lists what kind of boundary condition each boundary set is assigned to
  CArrayKokkos<int, array_layout, HostSpace, memory_traits> Boundary_Condition_Type_List;
  
  //number of displacement boundary conditions acting on nodes; used to size the reduced global stiffness map
  size_t Number_DOF_BCS;

  //! mapping used to get local ghost index from the global ID.
  //typedef ::Tpetra::Details::FixedHashTable<GO, LO, Kokkos::HostSpace::device_type>
    //global_to_local_table_host_type;

  //global_to_local_table_host_type global2local_map;
  //CArrayKokkos<int, Kokkos::LayoutLeft, Kokkos::HostSpace::device_type> active_ranks;

  //Pertains to local mesh information being stored as prescribed by the row map
  global_size_t local_nrows;
  global_size_t min_gid;
  global_size_t max_gid;
  global_size_t index_base;

  //allocation flags to avoid repeat MV and global matrix construction
  int Matrix_alloc;

  //debug flags
  int gradient_print_sync;

  //linear solver parameters
  Teuchos::RCP<Teuchos::ParameterList> Linear_Solve_Params;

  //multigrid solver data and functions
  //using equil_type = decltype (Tpetra::computeRowAndColumnOneNorms (*Global_Stiffness_Matrix, false));
  using equil_type = Tpetra::Details::EquilibrationInfo<real_t, device_type>;
  using mag_type = typename Kokkos::ArithTraits<real_t>::mag_type;
  using scaling_view_type = Kokkos::View<mag_type*, device_type>;
  equil_type equibResult;
  void equilibrateMatrix(Teuchos::RCP<Xpetra::Matrix<real_t,LO,GO,node_type> > &Axpetra, std::string equilibrate);
  void preScaleRightHandSides (Tpetra::MultiVector<real_t,LO,GO,node_type>& B, std::string equilibrate);
  void preScaleInitialGuesses (Tpetra::MultiVector<real_t,LO,GO,node_type>& X, std::string equilibrate);
  void postScaleSolutionVectors (Tpetra::MultiVector<real_t,LO,GO,node_type>& X, std::string equilibrate);
  void elementWiseMultiplyMultiVector (MV& X, scaling_view_type& scalingFactors,bool a);
  void elementWiseDivideMultiVector (MV& X, scaling_view_type& scalingFactors,bool a);
  void elementWiseMultiply (vec_array& X, scaling_view_type& scalingFactors, LO numRows, bool takeSquareRootsOfScalingFactors);
  void elementWiseDivide (vec_array& X, scaling_view_type& scalingFactors, LO numRows, bool takeSquareRootsOfScalingFactors);

  //division functor
  class ElementWiseDivide {
    public:
      static_assert (vec_array::Rank == 2, "ViewType1 must be a rank-2 "
                  "Kokkos::View in order to use this specialization.");

      ElementWiseDivide (const vec_array& X,
                         const scaling_view_type& scalingFactors, bool takeSquareRootsOfScalingflag) :
      X_ (X),
      scalingFactors_ (scalingFactors)
      { 
        takeAbsoluteValueOfScalingFactors = true;
        takeSquareRootsOfScalingFactors = takeSquareRootsOfScalingflag;
      }

      KOKKOS_INLINE_FUNCTION void operator () (const LO i) const {
        using val_type = typename scaling_view_type::non_const_value_type;
        using KAT = Kokkos::ArithTraits<val_type>;
        using mag_type = typename KAT::mag_type;
        using KAM = Kokkos::ArithTraits<mag_type>;

        for (LO j = 0; j < static_cast<LO> (X_.extent (1)); ++j) {
          if (takeAbsoluteValueOfScalingFactors) {
            const mag_type scalFactAbs = KAT::abs (scalingFactors_(i));
            const mag_type scalFinalVal = takeSquareRootsOfScalingFactors ?
              KAM::sqrt (scalFactAbs) : scalFactAbs;
            X_(i,j) = X_(i,j) / scalFinalVal;
          }
          else {
            const val_type scalFact = scalingFactors_(i);
            const val_type scalFinalVal = takeSquareRootsOfScalingFactors ?
              KAT::sqrt (scalFact) : scalFact;
            X_(i,j) = X_(i,j) / scalFinalVal;
          }
        }
      }

    private:
      vec_array X_;
      typename scaling_view_type::const_type scalingFactors_;
      bool takeAbsoluteValueOfScalingFactors;
      bool takeSquareRootsOfScalingFactors;
  };

  //multiplication functor
  class ElementWiseMultiply {
    public:
      static_assert (vec_array::Rank == 2, "ViewType1 must be a rank-2 "
                 "Kokkos::View in order to use this specialization.");

      ElementWiseMultiply (const vec_array& X,
                       const scaling_view_type& scalingFactors, bool takeSquareRootsOfScalingflag) :
      X_ (X),
      scalingFactors_ (scalingFactors)
      { 
        takeAbsoluteValueOfScalingFactors = true;
        takeSquareRootsOfScalingFactors = takeSquareRootsOfScalingflag;
      }

      KOKKOS_INLINE_FUNCTION void operator () (const LO i) const {
      using val_type = typename scaling_view_type::non_const_value_type;
      using KAT = Kokkos::ArithTraits<val_type>;
      using mag_type = typename KAT::mag_type;
      using KAM = Kokkos::ArithTraits<mag_type>;

      for (LO j = 0; j < static_cast<LO> (X_.extent (1)); ++j) {
        if (takeAbsoluteValueOfScalingFactors) {
          const mag_type scalFactAbs = KAT::abs (scalingFactors_(i));
          const mag_type scalFinalVal = takeSquareRootsOfScalingFactors ?
            KAM::sqrt (scalFactAbs) : scalFactAbs;
          X_(i,j) = X_(i,j) * scalFinalVal;
        }
        else {
          const val_type scalFact = scalingFactors_(i);
          const val_type scalFinalVal = takeSquareRootsOfScalingFactors ?
            KAT::sqrt (scalFact) : scalFact;
          X_(i,j) = X_(i,j) * scalFinalVal;
        }
      }
    }

    private:
      vec_array X_;
      typename scaling_view_type::const_type scalingFactors_;
      bool takeAbsoluteValueOfScalingFactors;
      bool takeSquareRootsOfScalingFactors;
  };
  
};

#endif // end HEADER_H
