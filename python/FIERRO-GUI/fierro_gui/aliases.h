#ifndef ALIASES_H
#define ALIASES_H
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

#include "host_types.h"
#include "kokkos_types.h"


using real_t = double;
using u_int  = unsigned int;


// simple data type names
namespace mtr
{
    // dense host (only cpu) types
    template <typename T>
    using CArrayHost = CArray <T>;

    template <typename T>
    using CMatrixHost = CMatrix <T>;

    template <typename T>
    using FArrayHost = FArray <T>;
    
    template <typename T>
    using FMatrixHost = FMatrix <T>;

    template <typename T>
    using ViewCArrayHost = ViewCArray <T>;
    
    template <typename T>
    using ViewCMatrixHost = ViewCMatrix <T>;
    
    template <typename T>
    using ViewFArrayHost = ViewFArray <T>;
    
    template <typename T>
    using ViewFMatrixHost = ViewFMatrix <T>;


    // ragged and sparse host (only cpu) types
    template <typename T>
    using RaggedCArrayHost = RaggedRightArray <T>;

    template <typename T>
    using RaggedFArrayHost = RaggedDownArray <T>;

    template <typename T>
    using DynamicRaggedCArrayHost = DynamicRaggedRightArray <T>;

    template <typename T>
    using DynamicRaggedFArrayHost = DynamicRaggedDownArray <T>;

    template <typename T>
    using CSRArrayHost = CSRArray <T>;

    template <typename T>
    using CSCArrayHost = CSCArray <T>;


} // end namespace


#ifdef HAVE_KOKKOS
namespace mtr
{
    // dense device types
    template <typename T>
    using CArrayDevice = CArrayKokkos <T>;

    template <typename T>
    using CMatrixDevice = CMatrixKokkos <T>;
    
    template <typename T>
    using FArrayDevice = FArrayKokkos <T>;
    
    template <typename T>
    using FMatrixDevice = FMatrixKokkos <T>;

    template <typename T>
    using ViewCArrayDevice = ViewCArrayKokkos <T>;
    
    template <typename T>
    using ViewCMatrixDevice = ViewCMatrixKokkos <T>;
    
    template <typename T>
    using ViewFArrayDevice = ViewFArrayKokkos <T>;
    
    template <typename T>
    using ViewFMatrixDevice = ViewFMatrixKokkos <T>;
    

    // ragged and sparse device types
    template <typename T>
    using RaggedCArrayDevice = RaggedRightArrayKokkos <T>;

    template <typename T>
    using RaggedFArrayDevice = RaggedDownArrayKokkos <T>;

    template <typename T>
    using DynamicRaggedCArrayDevice = DynamicRaggedRightArrayKokkos <T>;
    
    template <typename T>
    using DynamicRaggedFArrayDevice = DynamicRaggedDownArrayKokkos <T>;

    template <typename T>
    using CSRArrayDevice = CSRArrayKokkos <T>;
    
    template <typename T>
    using CSCArrayDevice = CSCArrayKokkos <T>;


    // dual dense types
    template <typename T>
    using CArrayDual = DCArrayKokkos <T>;
    
    template <typename T>
    using CMatrixDual = DCMatrixKokkos <T>;
    
    template <typename T>
    using FArrayDual = DFArrayKokkos <T>;
    
    template <typename T>
    using FMatrixDual = DFMatrixKokkos <T>;

    template <typename T>
    using ViewCArrayDual = DViewCArrayKokkos <T>;
    
    template <typename T>
    using ViewCMatrixDual = DViewCMatrixKokkos <T>;
    
    template <typename T>
    using ViewFArrayDual = DViewFArrayKokkos <T>;
    
    template <typename T>
    using ViewFMatrixDual = DViewFMatrixKokkos <T>;


} // end namespace
#endif // end if have Kokkos for simple data type names

#endif // ALIASES_H
