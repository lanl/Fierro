#include "evpfft.h"
#include "XdmfUniformGridWriter.h"
#include "math_functions.h"
#include <stdio.h>
#include <array>
#include "Profiler.h"
#include "vm.h"

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

        sprintf (dataset_name, "%s_%d%d", data_name, ii, jj);        
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

template <typename T>
void write_N1xN2xN3(T *data_ptr, int npts1, int npts2, int npts3,
  const char *data_name, const char* group_name, std::shared_ptr<MicroOutputWriter<MPI_ORDER_FORTRAN>> micro_writer,
  XdmfUniformGridWriter &xdmf_writer, const char* hdf5_filename, int my_rank)
{
  /* This function assumes data_ptr has a Fortran Array Order */

  ViewFMatrix <T> data_view (data_ptr,npts1,npts2,npts3);
  char dataset_full_path [200];
  std::string dataset_name = data_name;

  // for xdmf_writer
  char aLocation [200];
  std::string aNumberType;
  std::string aPrecision;
  xdmf_writer.deduce_attribute_number_type <T> (data_ptr, aNumberType, aPrecision);

  sprintf (dataset_full_path, "%s/%s", group_name, dataset_name.c_str());
  micro_writer->write <T> (dataset_full_path, data_ptr);

  if (0 == my_rank) {
    sprintf (aLocation, "%s:%s", hdf5_filename, dataset_full_path);
    xdmf_writer.write_attribute(dataset_name.c_str(), "Scalar", aNumberType.c_str(), aPrecision.c_str(), aLocation);
  }

}

void calculate_vm_field(MatrixTypeRealDual field, MatrixTypeRealDual vm_field, int vm_type)
{
 /* 
 * vm_type 0: stress 
 * vm_type 1: strain
 * */

  if (vm_type != 0 && vm_type != 1) {
    throw std::invalid_argument("vm_type must be 0 or 1");
  }

  int npts1 = field.dims(3);
  int npts2 = field.dims(4);
  int npts3 = field.dims(5);

  FOR_ALL(
          k, 1, npts3+1,
          j, 1, npts2+1,
          i, 1, npts1+1, {

    real_t aux_[3*3];
    ViewMatrixTypeReal aux(aux_,3,3);

    for (int ii = 1; ii <= 3; ii++) {
      for (int jj = 1; jj <= 3; jj++) {
        aux(ii,jj) = field(ii,jj,i,j,k);
      }
    }
    
    if (vm_type == 0) {
      vm_field(i,j,k) = vm_stress(aux.pointer());
    }
    else if (vm_type == 1) {
      vm_field(i,j,k) = vm(aux.pointer());
    }

  }); // end FOR_ALL
  Kokkos::fence();

  vm_field.update_host();
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

        invert_matrix <6> (cg66aux.pointer());
          
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

  eel.update_host();
}

void symmetrize(
  const MatrixTypeRealDual& unsymmetrized_field,
  MatrixTypeRealDual& symmetrized_field)
{
  size_t npts1 = unsymmetrized_field.dims(3);
  size_t npts2 = unsymmetrized_field.dims(4);
  size_t npts3 = unsymmetrized_field.dims(5);

  FOR_ALL(
          k, 1, npts3+1,
          j, 1, npts2+1,
          i, 1, npts1+1, {

       for (int ii = 1; ii <= 3; ii++) {
          for (int jj = 1; jj <= 3; jj++) {
            symmetrized_field(ii,jj,i,j,k) = 
              0.5*(unsymmetrized_field(ii,jj,i,j,k) + unsymmetrized_field(jj,ii,i,j,k));
          }   
        }

  }); // end FOR_ALL
  Kokkos::fence();

  symmetrized_field.update_host(); 
  return; 
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

  // calculate elastic-strain
  MatrixTypeRealDual eel (3,3,npts1,npts2,npts3);
  calculate_eel(eel);

  // Some arrays for results
  MatrixTypeRealDual sym_disgrad (3,3,npts1,npts2,npts3);
  MatrixTypeRealDual sym_edotp (3,3,npts1,npts2,npts3);
  MatrixTypeRealDual sym_eel (3,3,npts1,npts2,npts3);

  // Symmetrize some fields
  symmetrize(disgrad, sym_disgrad);
  symmetrize(edotp, sym_edotp);
  symmetrize(eel, sym_eel);
  
  // write strain field
  disgrad.update_host();
  write_distinct_components_of_3x3xN1xN2xN3 <real_t> (sym_disgrad.host_pointer(), npts1, npts2, npts3, "strain",
                               group_name, micro_writer, xdmf_writer, micro_writer->filename_, my_rank);

  // write stress field
  sg.update_host();
  write_distinct_components_of_3x3xN1xN2xN3 <real_t> (sg.host_pointer(), npts1, npts2, npts3, "stress",
                               group_name, micro_writer, xdmf_writer, micro_writer->filename_, my_rank);

  // write plastic-strain-rate field
  edotp.update_host();
  write_distinct_components_of_3x3xN1xN2xN3 <real_t> (sym_edotp.host_pointer(), npts1, npts2, npts3, "pl-strain-rate",
                               group_name, micro_writer, xdmf_writer, micro_writer->filename_, my_rank);

  // write elastic-strain
  write_distinct_components_of_3x3xN1xN2xN3 <real_t> (sym_eel.host_pointer(), npts1, npts2, npts3, "el-strain",
                               group_name, micro_writer, xdmf_writer, micro_writer->filename_, my_rank);

  // write jphase
  write_N1xN2xN3 <int> (jphase.host_pointer(), npts1, npts2, npts3, "phase_id",
                        group_name, micro_writer, xdmf_writer, micro_writer->filename_, my_rank);

  // for writing vm fields 
  MatrixTypeRealDual vm_field (npts1,npts2,npts3);

  // calculate and write strain vm field
  calculate_vm_field(sym_disgrad, vm_field, 1);
  write_N1xN2xN3 <real_t> (vm_field.host_pointer(), npts1, npts2, npts3, "vm-strain",
                        group_name, micro_writer, xdmf_writer, micro_writer->filename_, my_rank);

  // calculate and write stress vm field
  calculate_vm_field(sg, vm_field, 0);
  write_N1xN2xN3 <real_t> (vm_field.host_pointer(), npts1, npts2, npts3, "vm-stress",
                        group_name, micro_writer, xdmf_writer, micro_writer->filename_, my_rank);

  // calculate and write plastic-strain-rate vm field
  calculate_vm_field(sym_edotp, vm_field, 1);
  write_N1xN2xN3 <real_t> (vm_field.host_pointer(), npts1, npts2, npts3, "vm-pl-strain-rate",
                        group_name, micro_writer, xdmf_writer, micro_writer->filename_, my_rank);

  // calculate and write elastic-strain vm field
  calculate_vm_field(sym_eel, vm_field, 1);
  write_N1xN2xN3 <real_t> (vm_field.host_pointer(), npts1, npts2, npts3, "vm-el-strain",
                        group_name, micro_writer, xdmf_writer, micro_writer->filename_, my_rank);



  // close xdmf file
  if (0 == my_rank) {
    xdmf_writer.write_footer();
    xdmf_writer.close_file();
  }
    
  // close group
  herr_t status = H5Gclose(group_id);
}
