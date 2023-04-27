#include "evpfft.h"
#include "VTKUniformGridWriter.h"
#include "inverse.h"
#include <stdio.h>
#include <array>

template <typename T>
void write_distinct_components_of_3x3xN1xN2xN3(T *data_ptr,
                                               int npts1, int npts2, int npts3,
                                               const char *data_name,
                                               VTKUniformGridWriter &vtk_writer);


template <typename T>
void write_distinct_components_of_3x3xN1xN2xN3(T *data_ptr,
                                               int npts1, int npts2, int npts3,
                                               const char *data_name,
                                               VTKUniformGridWriter &vtk_writer)
{
  /* This function assumes data_ptr has a Fortran Array Order */

  ViewFMatrix <T> data_view (data_ptr,3,3,npts1,npts2,npts3);
  FMatrix <T> buffer (npts1,npts2,npts3);
  char dataset_name [50];

  for (int ii = 1; ii <= 3; ii++) {
    for (int jj = 1; jj <= 3; jj++) {

      if (ii <= jj) {
        
        for (int k = 1; k <= npts3; k++) {
          for (int j = 1; j <= npts2; j++) {
            for (int i = 1; i <= npts1; i++) {         
              buffer(i,j,k) = data_view(ii,jj,i,j,k);
            } // end for k
          } // end for j
        } // end for i

        sprintf (dataset_name, "%s%d%d", data_name, ii, jj);        
        vtk_writer.write_scalar_dataset <T> (buffer.pointer(), dataset_name, "double");

      } // end if (ii <= jj)

    } // end for jj
  } // end for ii
}

void EVPFFT::calculate_eel(MatrixTypeRealDual &eel)
{
  FOR_ALL_CLASS(k, 1, npts3+1,
          j, 1, npts2+1,
          i, 1, npts1+1, {

      int jph;
      real_t cg66aux_ [6*6];
      real_t cinvg_ [3*3*3*3];
      ViewMatrixTypeReal cg66aux(cg66aux_,6,6);
      ViewMatrixTypeReal cinvg(cinvg_,3,3,3,3);

      jph = jphase(i,j,k);

      if (igas(jph)  == 0) {

        for (int jj = 1; jj <= 6; jj++) {
          for (int ii = 1; ii <= 6; ii++) {
            cg66aux(ii,jj) = cg66(ii,jj,i,j,k);
          }
        }

#ifdef LU_MATRIX_INVERSE
        lu_inverse(cg66aux.pointer(), 6);
#elif GJE_MATRIX_INVERSE
        inverse_gj(cg66aux.pointer(), 6);
#endif
          
        cb.chg_basis_3(cg66aux.pointer(), cinvg.pointer(), 3, 6, cb.B_basis_device_pointer());

        for (int jj = 1; jj <= 3; jj++) {
          for (int ii = 1; ii <= 3; ii++) {
            eel(ii,jj,i,j,k) = 0.0;
          }
        }

        for (int kk = 1; kk <= 3; kk++) {
          for (int ll = 1; ll <= 3; ll++) {
            for (int jj = 1; jj <= 3; jj++) {
              for (int ii = 1; ii <= 3; ii++) {
                eel(ii,jj,i,j,k) += cinvg(ii,jj,kk,ll) * sg(kk,ll,i,j,k);
              }
            }
          }
        }

      } else { // igas else

        for (int jj = 1; jj <= 3; jj++) {
          for (int ii = 1; ii <= 3; ii++) {
            eel(ii,jj,i,j,k) = -100.0;
          }
        }

      } // end if (igas(jph)  == 0)

  }); // end FOR_ALL_CLASS
  Kokkos::fence();
}

void EVPFFT::write_micro_state(int imicro)
{
if (iwfields == 1 and imicro % iwstep == 0) {

  // vtk file stuff
  VTKUniformGridWriter vtk_writer(npts3, npts2, npts1);
  char buffer [200];
  sprintf (buffer, "micro_state_timestep_%d.vtk", imicro); 
  vtk_writer.open_file(buffer);
  vtk_writer.write_header();
 
  // write strain field
  disgrad.update_host();
  write_distinct_components_of_3x3xN1xN2xN3 <real_t> (disgrad.host_pointer(),
                                                      npts1, npts2, npts3, "e",
                                                      vtk_writer);

  // write stress field
  sg.update_host();
  write_distinct_components_of_3x3xN1xN2xN3 <real_t> (sg.host_pointer(),
                                                      npts1, npts2, npts3, "s",
                                                      vtk_writer);

  // write plastic-strain-rate field
  edotp.update_host();
  write_distinct_components_of_3x3xN1xN2xN3 <real_t> (edotp.host_pointer(), 
                                                      npts1, npts2, npts3, "edotp",
                                                      vtk_writer);

  // calculate and write elastic-strain
  MatrixTypeRealDual eel (3,3,npts1,npts2,npts3);
  calculate_eel(eel);
  eel.update_host();
  write_distinct_components_of_3x3xN1xN2xN3 <real_t> (eel.host_pointer(), 
                                                      npts1, npts2, npts3, "eel",
                                                      vtk_writer);
  // close vtk file
  vtk_writer.close_file();
    
} // end if (iwfields and imicro % iwstep == 0)
}
