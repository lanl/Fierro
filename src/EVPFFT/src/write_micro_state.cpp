#include "evpfft.h"
#include "XdmfUniformGridWriter.h"
#include "math_functions.h"
#include <stdio.h>
#include <array>
#include "Profiler.h"

template <typename T>
void write_distinct_components_of_3x3xN1xN2xN3(T *data_ptr, int npts1, int npts2, int npts3,
  const char *data_name, const char* group_name, std::shared_ptr<MicroOutputWriter<MPI_ORDER_FORTRAN>> micro_writer,
  XdmfUniformGridWriter &xdmf_writer, const char* hdf5_filename, int my_rank);


template <typename T>
void write_distinct_components_of_3x3xN1xN2xN3(T *data_ptr, int npts1, int npts2, int npts3,
  const char *data_name, const char* group_name, std::shared_ptr<MicroOutputWriter<MPI_ORDER_FORTRAN>> micro_writer,
  XdmfUniformGridWriter &xdmf_writer, const char* hdf5_filename, int my_rank)
{
  /* This function assumes data_ptr has a Fortran Array Order */

  ViewFMatrix <T> data_view (data_ptr,3,3,npts1,npts2,npts3);
  FMatrix <T> buffer (npts1,npts2,npts3);
  char dataset_name [50];
  char dataset_full_path [200];

  // for xdmf_writer
  char aLocation [200];
  std::string aNumberType;
  std::string aPrecision;
  xdmf_writer.deduce_attribute_number_type <T> (data_ptr, aNumberType, aPrecision);

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
        sprintf (dataset_full_path, "%s/%s", group_name, dataset_name);
        micro_writer->write <T> (dataset_full_path, buffer.pointer());

        if (0 == my_rank) {
          sprintf (aLocation, "%s:%s", hdf5_filename, dataset_full_path);
          xdmf_writer.write_attribute(dataset_name, "Scalar", aNumberType.c_str(), aPrecision.c_str(), aLocation);
        }

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

        invert_matrix(cg66aux.pointer(), 6);
          
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

void EVPFFT::write_micro_state()
{
  Profiler profiler(__FUNCTION__);

  // xdmf file stuff
  XdmfUniformGridWriter xdmf_writer(npts3_g, npts2_g, npts1_g);
  if (0 == my_rank) {
    char buffer [100];

    sprintf (buffer, "micro_state_timestep_%d.xdmf", imicro); 
    xdmf_writer.open_file(buffer);

    sprintf (buffer, "TimeStep_%d", imicro);
    xdmf_writer.write_header(buffer);
  }
 
  // create group
  char group_name [50];
  sprintf (group_name, "/TimeStep_%d", imicro);  
  hid_t group_id = H5Gcreate(micro_writer->file_id_, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // write strain field
  disgrad.update_host();
  write_distinct_components_of_3x3xN1xN2xN3 <real_t> (disgrad.host_pointer(), npts1, npts2, npts3, "e",
                               group_name, micro_writer, xdmf_writer, micro_writer->filename_, my_rank);

  // write stress field
  sg.update_host();
  write_distinct_components_of_3x3xN1xN2xN3 <real_t> (sg.host_pointer(), npts1, npts2, npts3, "s",
                               group_name, micro_writer, xdmf_writer, micro_writer->filename_, my_rank);

  // write plastic-strain-rate field
  edotp.update_host();
  write_distinct_components_of_3x3xN1xN2xN3 <real_t> (edotp.host_pointer(), npts1, npts2, npts3, "edotp",
                               group_name, micro_writer, xdmf_writer, micro_writer->filename_, my_rank);

  // calculate and write elastic-strain
  MatrixTypeRealDual eel (3,3,npts1,npts2,npts3);
  calculate_eel(eel);
  eel.update_host();
  write_distinct_components_of_3x3xN1xN2xN3 <real_t> (eel.host_pointer(), npts1, npts2, npts3, "eel",
                               group_name, micro_writer, xdmf_writer, micro_writer->filename_, my_rank);
  // close xdmf file
  if (0 == my_rank) {
    xdmf_writer.write_footer();
    xdmf_writer.close_file();
  }
    
  // close group
  herr_t status = H5Gclose(group_id);
}
