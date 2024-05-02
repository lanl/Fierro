// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef FIERRO_MMA_OBJECTIVE_H
#define FIERRO_MMA_OBJECTIVE_H

#include <string>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include "Tpetra_Details_DefaultTypes.hpp"

#include "ROL_Types.hpp"
#include <ROL_TpetraMultiVector.hpp>
#include "ROL_Elementwise_Reduce.hpp"
#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"

/** @ingroup func_group
    \class ROL::ObjectiveMMA
    \brief Provides the interface to to Method of Moving Asymptotes 
           Objective function

    ---
*/

class ObjectiveMMA : public ROL::Objective<real_t> {
  
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Tpetra::Map<>::node_type Node;
  typedef Tpetra::Map<LO, GO, Node> Map;
  typedef Tpetra::MultiVector<real_t, LO, GO, Node> MV;
  typedef ROL::Vector<real_t> V;
  typedef const ROL::Vector<real_t> const_V;
  typedef ROL::TpetraMultiVector<real_t,LO,GO,Node> ROL_MV;
  
  using traits = Kokkos::ViewTraits<LO*, Kokkos::LayoutLeft, void, void>;
  using array_layout    = typename traits::array_layout;
  using execution_space = typename traits::execution_space;
  using device_type     = typename traits::device_type;
  using memory_traits   = typename traits::memory_traits;
  using global_size_t = Tpetra::global_size_t;

  typedef Kokkos::View<real_t*, Kokkos::LayoutRight, device_type, memory_traits> values_array;
  typedef Kokkos::View<GO*, array_layout, device_type, memory_traits> global_indices_array;
  typedef Kokkos::View<LO*, array_layout, device_type, memory_traits> indices_array;
  
  //typedef Kokkos::DualView<real_t**, Kokkos::LayoutLeft, device_type>::t_dev vec_array;
  typedef MV::dual_view_type::t_dev vec_array;
  typedef MV::dual_view_type::t_host host_vec_array;
  typedef Kokkos::View<const real_t**, array_layout, HostSpace, memory_traits> const_host_vec_array;
  typedef MV::dual_view_type dual_vec_array;
  typedef ROL::Objective<real_t>       OBJ;
  typedef ROL::BoundConstraint<real_t> BND;
  int update_count = 0;

private:

  const ROL::Ptr<OBJ> obj_;
  const ROL::Ptr<BND> bnd_;
  

  ROL::Ptr<const_V> l_shift; // Lower bound
  ROL::Ptr<const_V> u_shift; // Upper bound
  ROL::Ptr<V> lv_; // Lower bound
  ROL::Ptr<V> uv_; // Upper bound
  
  ROL::Ptr<V> p_; // First MMA numerator
  ROL::Ptr<V> q_; // Second MMA numerator

  ROL::Ptr<V> d_; // Scratch vector

  real_t fval_; // Original objective value

  real_t tol_;

   ROL::Ptr<const MV> getVector( const V& x ) {
    return dynamic_cast<const ROL_MV&>(x).getVector();
  }

  ROL::Ptr<MV> getVector( V& x ) {
    return dynamic_cast<ROL_MV&>(x).getVector();
  }

public:

  ObjectiveMMA( const ROL::Ptr<ROL::Objective<real_t> > &obj,
                const ROL::Ptr<ROL::BoundConstraint<real_t> > &bnd,
                const ROL::Ptr<V> &x,
                real_t tol=std::sqrt(ROL::ROL_EPSILON<real_t>()) ) : 
    obj_(obj), bnd_(bnd), tol_(tol) {

    l_shift   = bnd_->getLowerBound();
    u_shift   = bnd_->getUpperBound();
    

    p_   = x->clone();
    q_   = x->clone();
    d_   = x->clone();
    lv_  = x->clone();
    uv_  = x->clone();

    update_count = 0;
 
  }

  /* --------------------------------------------------------------------------------------
   Update solver state variables to synchronize with the current design variable vector, z
  ----------------------------------------------------------------------------------------- */

  void update(const ROL::Vector<real_t> &x, ROL::UpdateType type, int iter = -1 ) {
  
    ROL::Elementwise::ThresholdUpper<real_t> positive(0.0);
    ROL::Elementwise::ThresholdLower<real_t> negative(0.0);
    ROL::Elementwise::Power<real_t>          square(2.0);
    ROL::Elementwise::Multiply<real_t>       mult;
    ROL::Elementwise::Scale<real_t>          negative_one(-1);
    ROL::Elementwise::Scale<real_t>          limit_scalar(0.25);
    ROL::Elementwise::Scale<real_t>          limit_invscalar(4);

    if (type == ROL::UpdateType::Initial)  {
      // This is the first call to update
      obj_->update(x,type,iter);
      
      fval_ = obj_->value(x,tol_);
      obj_->gradient(*p_,x,tol_);
      q_->set(*p_);
      lv_->set(x);
      uv_->set(x);
      lv_->axpy(1.0,*l_shift);
      uv_->axpy(1.0,*u_shift);

      p_->applyUnary(positive);
      q_->applyUnary(negative);

      d_->set(*uv_);
      d_->axpy(-1.0,x);
      d_->applyUnary(square);
      p_->applyBinary(mult,*d_);

      d_->set(x);
      d_->axpy(-1.0,*lv_);
      d_->applyUnary(square);
      q_->applyBinary(mult,*d_);
      q_->applyUnary(negative_one);

    }
    else if (type == ROL::UpdateType::Accept) {
      // u_ was set to u=S(x) during a trial update
      // and has been accepted as the new iterate
      // This is the first call to update
      obj_->update(x,type,iter);
      
      fval_ = obj_->value(x,tol_);
      obj_->gradient(*p_,x,tol_);
      q_->set(*p_);
      lv_->set(x);
      uv_->set(x);
      lv_->axpy(1.0,*l_shift);
      uv_->axpy(1.0,*u_shift);

      p_->applyUnary(positive);
      q_->applyUnary(negative);

      d_->set(*uv_);
      d_->axpy(-1.0,x);
      d_->applyUnary(square);
      p_->applyBinary(mult,*d_);

      d_->set(x);
      d_->axpy(-1.0,*lv_);
      d_->applyUnary(square);
      q_->applyBinary(mult,*d_);
      q_->applyUnary(negative_one);
    }
    else if (type == ROL::UpdateType::Revert) {
      // u_ was set to u=S(x) during a trial update
      // and has been rejected as the new iterate
      // Revert to cached value
      // This is the first call to update
      // u_ was set to u=S(x) during a trial update
      // and has been accepted as the new iterate
      // This is the first call to update
      obj_->update(x,type,iter);
      
      fval_ = obj_->value(x,tol_);
      obj_->gradient(*p_,x,tol_);
      q_->set(*p_);
      lv_->set(x);
      uv_->set(x);
      lv_->axpy(1.0,*l_shift);
      uv_->axpy(1.0,*u_shift);

      p_->applyUnary(positive);
      q_->applyUnary(negative);

      d_->set(*uv_);
      d_->axpy(-1.0,x);
      d_->applyUnary(square);
      p_->applyBinary(mult,*d_);

      d_->set(x);
      d_->axpy(-1.0,*lv_);
      d_->applyUnary(square);
      q_->applyBinary(mult,*d_);
      q_->applyUnary(negative_one);
    }
    else if (type == ROL::UpdateType::Trial) {
      // This is a new value of x
      obj_->update(x,type,iter);
      
      fval_ = obj_->value(x,tol_);
      obj_->gradient(*p_,x,tol_);
      q_->set(*p_);
      lv_->set(x);
      uv_->set(x);
      lv_->axpy(1.0,*l_shift);
      uv_->axpy(1.0,*u_shift);

      p_->applyUnary(positive);
      q_->applyUnary(negative);

      d_->set(*uv_);
      d_->axpy(-1.0,x);
      d_->applyUnary(square);
      p_->applyBinary(mult,*d_);

      d_->set(x);
      d_->axpy(-1.0,*lv_);
      d_->applyUnary(square);
      q_->applyBinary(mult,*d_);
      q_->applyUnary(negative_one);

      // This is the first call to update
    }
    else { // ROL::UpdateType::Temp
      // This is a new value of x used for,
      // e.g., finite-difference checks
      
      // This is the first call to update
    }

  }
  
  /* ------------------------------------------------------------------------------------------------------
   Update objective value with the current design variable vector, z
  \f[ F(x) \approx F(x^0) + \sum\limit_{i=1}^n \left( \frac{p_i}{U_i-x_i} + \frac{q_i}{x_i-L_i}\right) \f] 
  --------------------------------------------------------------------------------------------------------- */

  real_t value( const ROL::Vector<real_t> &x, real_t &tol ) {
  
    ROL::Elementwise::ReductionSum<real_t>    sum;
    ROL::Elementwise::DivideAndInvert<real_t> divinv;
    real_t fval = fval_;
    
    //debug print of design variables
    //ROL::Ptr<MV> pp = getVector(*p_);
    //std::ostream &out = std::cout;
    //Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    //(*fos).setOutputToRootOnly(0);
    //*fos << "Gradient data :" << std::endl;
    //pp->describe(*fos,Teuchos::VERB_EXTREME);

    d_->set(*uv_);
    d_->axpy(-1.0,x);
    d_->applyBinary(divinv,*p_);  

    fval += d_->reduce(sum);
    d_->set(x);
    d_->axpy(-1.0,*lv_);
    d_->applyBinary(divinv,*q_);

    fval += d_->reduce(sum);   
    
    return fval;

  }

  /* --------------------------------------------------------------------------------------
   Update gradient vector (g) with the current design variable vector, z
     \f[ \frac{F(x)}{\partial x_j} =  \frac{p_j}{(U_j-x_j)^2} - \frac{q_j}({x_j-L_j)^2}\ \f]
  ----------------------------------------------------------------------------------------- */

  void gradient( ROL::Vector<real_t> &g, const ROL::Vector<real_t> &x, real_t &tol ) {

    ROL::Elementwise::DivideAndInvert<real_t> divinv;
    ROL::Elementwise::Power<real_t>           square(2.0);
    ROL::Elementwise::Scale<real_t>          negative_one(-1);

    d_->set(*uv_);
    d_->axpy(-1.0,x);
    d_->applyUnary(square);
    d_->applyBinary(divinv,*p_);            

    g.set(*d_);

    d_->set(x);
    d_->axpy(-1.0,*lv_);
    d_->applyUnary(square);
    d_->applyBinary(divinv,*q_);
    d_->applyUnary(negative_one);

    g.plus(*d_);

  }

  /* --------------------------------------------------------------------------------------
   Update Hessian vector product (hv) using the differential design vector (v) and
   the current design variable vector, z
  ----------------------------------------------------------------------------------------- */
  
  void hessVec( ROL::Vector<real_t> &hv, const ROL::Vector<real_t> &v, const ROL::Vector<real_t> &x, real_t &tol ) {
  
    ROL::Elementwise::DivideAndInvert<real_t> divinv;
    ROL::Elementwise::Multiply<real_t>        mult;
    ROL::Elementwise::Power<real_t>           cube(3.0);

    d_->set(*uv_);
    d_->axpy(-1.0,x);
    d_->applyUnary(cube);          
    d_->applyBinary(divinv,*p_);           
    d_->scale(2.0);

    hv.set(*d_); 

    d_->set(x);
    d_->axpy(-1.0,*lv_);
    d_->applyUnary(cube);
    d_->applyBinary(divinv,*q_);
    d_->scale(2.0);
    
    hv.plus(*d_);
    hv.applyBinary(mult,v);

  }

  void invHessVec( ROL::Vector<real_t> &hv, const ROL::Vector<real_t> &v, const ROL::Vector<real_t> &x, real_t &tol ) {

    ROL::Elementwise::DivideAndInvert<real_t> divinv;
    ROL::Elementwise::Multiply<real_t>        mult;
    ROL::Elementwise::Power<real_t>           cube(3.0);

    d_->set(*uv_);
    d_->axpy(-1.0,x);
    d_->applyUnary(cube);          
    d_->applyBinary(divinv,*p_);           
    d_->scale(2.0);

    hv.set(*d_); 

    d_->set(x);
    d_->axpy(-1.0,*lv_);
    d_->applyUnary(cube);
    d_->applyBinary(divinv,*q_);
    d_->scale(2.0);
    
    hv.plus(*d_);
    hv.applyBinary(divinv,v);
      
  }

}; // class ObjectiveMMA 


#endif // FIERRO_MMA_OBJECTIVE

