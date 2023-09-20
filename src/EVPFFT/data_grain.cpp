#include <stdio.h>
#include "evpfft.h"
#include "definitions.h"
#include "utilities.h"
#include "euler.h"
#include "inverse.h"
#include "voigt.h"

using namespace utils;


void EVPFFT::data_grain(const std::string & filetext)
{
  int nph1 = 0;

  switch (cmd.micro_filetype)
  {
    case 0:
    read_classic_los_alamos_texture_file(filetext, nph1);
    break;

    case 1:
    read_hdf5_texture_file(filetext, nph1);
    break;

    default:
      throw std::runtime_error("invalid micro_filetype.");
  }

  for (int kk = 1; kk <= npts3; kk++) {
    for (int jj = 1; jj <= npts2; jj++) {
      for (int ii = 1; ii <= npts1; ii++) {
                    
        int jph;
        real_t ph, th, om;
        real_t ph_rad, th_rad, om_rad;
        real_t dum;
      
        real_t aa_[3*3];
        real_t caux3333_[3*3*3*3];
        real_t caux66_[6*6];
        ViewMatrixTypeReal aa(aa_,3,3);
        ViewMatrixTypeReal caux3333(caux3333_,3,3,3,3);
        ViewMatrixTypeReal caux66(caux66_,6,6);
 
        ph  = ph_array(ii,jj,kk);
        th  = th_array(ii,jj,kk);
        om  = om_array(ii,jj,kk);
        jph = jphase.host(ii,jj,kk);
           
        // CALCULATES THE TRANSFORMATION MATRIX AA WHICH TRANSFORMS FROM
        // SAMPLE TO CRYSTAL. STORES AG, WHICH TRANSFORMS FROM CRYSTAL TO SAMPLE.
        ph_rad = ph*pi/180.0;
        th_rad = th*pi/180.0;
        om_rad = om*pi/180.0;
        euler(2, ph_rad, th_rad, om_rad, aa.pointer());
                    
        for (int j = 1; j <= 3; j++) {
          for (int k = 1; k <= 3; k++) {
            ag.host(j,k,ii,jj,kk) = aa(k,j);
          }
        }

        for (int i1 = 1; i1 <= 3; i1++) {
          for (int j1 = 1; j1 <= 3; j1++) {
            for (int k1 = 1; k1 <= 3; k1++) {
              for (int l1 = 1; l1 <= 3; l1++) {
            
                dum = 0.0;
              
                for (int i2 = 1; i2 <= 3; i2++) {
                  for (int j2 = 1; j2 <= 3; j2++) {
                    for (int k2 = 1; k2 <= 3; k2++) {
                      for (int l2 = 1; l2 <= 3; l2++) {
                    
                        dum += aa(i2,i1) *
                               aa(j2,j1) *
                               aa(k2,k1) *
                               aa(l2,l1) *
                               cc(i2,j2,k2,l2,jph);
                             
                      }
                    }
                  }
                }
                            
                caux3333(i1,j1,k1,l1) = dum;
             
              }
            }
          }
        }
            
      cb.chg_basis_4(caux66.pointer(), caux3333.pointer(), 4, 6, cb.B_basis_host_pointer());
      
      for (int i = 1; i <= 6; i++) {
        for (int j = 1; j <= 6; j++) {
          cg66.host(i,j,ii,jj,kk) = caux66(i,j);
          
          // FIXME, this is racy, needs proper reduction (Original Fortran code comment)
          //Kokkos::atomic_add(&c066(i,j), caux66(i,j)*wgt);
          c066.host(i,j) += caux66(i,j)*wgt;
        }
      }
     
      } // end for ii
    } // end for jj
  } // end for kk

  // update device
  c066.update_device();
  cg66.update_device();
  ag.update_device();
  Kokkos::fence();

#if 0
  // debug print
  print_array_to_file(c066.host_pointer(), c066.size(), my_rank, "c066.txt");
#endif

  wph1 = (1.0*nph1)/(npts1*npts2*npts3);

  MatrixTypeRealHost s066(6,6);
  for (int i = 1; i <= 6; i++) {
    for (int j = 1; j <= 6; j++) {
      s066(i,j) = c066.host(i,j);
    }
  }

#ifdef LU_MATRIX_INVERSE
  lu_inverse(s066.pointer(),6);
#elif GJE_MATRIX_INVERSE
  inverse_gj(s066.pointer(),6);
#endif

  cb.chg_basis_3(c066.host_pointer(), c0.host_pointer(), 3, 6, cb.B_basis_host_pointer());
  cb.chg_basis_3(s066.pointer(), s0.host_pointer(), 3, 6, cb.B_basis_host_pointer());

  // update device
  c0.update_device();
  s0.update_device();
  Kokkos::fence();
}


void EVPFFT::read_classic_los_alamos_texture_file(const std::string & filetext, int & nph1)
{
  std::ifstream ur2;
  ur2.open(filetext);
  check_that_file_is_open(ur2, filetext.c_str());

  real_t ph, th, om;
  int iiii, jjjj, kkkk, jgr, jph;

  for (int kk = 1; kk <= npts3; kk++) {
    for (int jj = 1; jj <= npts3; jj++) {
      for (int ii = 1; ii <= npts3; ii++) {
          
        ur2 >> ph >> th >> om >> iiii >> jjjj >> kkkk >> jgr >> jph; CLEAR_LINE(ur2);
            
        ph_array(ii,jj,kk) = ph; 
        th_array(ii,jj,kk) = th; 
        om_array(ii,jj,kk) = om; 
        jgrain(ii,jj,kk)   = jgr;
        jphase.host(ii,jj,kk) = jph;
            
        if (jph == 1) nph1 = nph1 + 1;
        if (igas.host(jph) == 1) exit(1); // this got broken when the read loop got split
          
      }   
    }   
  }
  // update device
  jphase.update_device();
  Kokkos::fence();
  ur2.close();
}

void EVPFFT::read_hdf5_texture_file(const std::string & filetext, int & nph1)
{
  throw std::runtime_error("EVPFFT::read_hdf5_texture_file(...) not yet implemented");
}

