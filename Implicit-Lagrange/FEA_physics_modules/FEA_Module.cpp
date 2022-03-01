#include "FEA_Module.h"
using namespace utils;

FEA_Module::FEA_Module(){
  num_nodes = 0;
  hessvec_count = update_count = 0;
  file_index = 0;
  linear_solve_time = hessvec_time = hessvec_linear_time = 0;

  Matrix_alloc=0;
  gradient_print_sync = 0;
  //RCP pointer to *this (Parallel Nonlinear Solver Object)
  //FEM_pass = Teuchos::rcp(this);

  //preconditioner construction
  Hierarchy_Constructed = false;

  //Trilinos output stream
  std::ostream &out = std::cout;
  fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
  (*fos).setOutputToRootOnly(0);
}

FEA_Module::~FEA_Module() {}