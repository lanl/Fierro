#ifndef MATAR_H
#define MATAR_H
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


//  ----
//   Standard data structures (host types)
//   1. FArray
//   2. ViewFArray
//   3. FMatrix
//   4. ViewFMatrix
//   5. CArray
//   6. ViewCArray
//   7. CMatrix
//   8. ViewCMatrix
//   9. RaggedRightArray
//   10. RaggedDownArray
//   11. DynamicRaggedRightArray
//   12. DynamicRaggedDownArray
//   13. CSRArray
//   14. CSCArray

//  ----
//   Kokkos Data structures (device and dual types)
//   15. FArrayKokkos
//   16. ViewFArrayKokkos
//   17. FMatrixKokkos
//   18. ViewFMatrixKokkos
//   19. CArrayKokkos
//   20. ViewCArrayKokkos
//   21. CMatrixKokkos
//   22. ViewCMatrixKokkos
//   23. RaggedRightArrayKokkos
//   24. RaggedDownArrayKokkos
//   25. DynamicRaggedRightArrayKokkos
//   26. DynamicRaggedDownArrayKokkos
//   27. CSRArrayKokkos
//   28. CSCArrayKokkos
//   29. DViewCArrayKokkos
//   30. DViewCMatrixKokkos
//   31. DViewFArrayKokkos
//   32. DViewFMatrixKokkos


#include "macros.h"
#include "host_types.h"
#include "kokkos_types.h"
#include "aliases.h"



#endif // MATAR_H
