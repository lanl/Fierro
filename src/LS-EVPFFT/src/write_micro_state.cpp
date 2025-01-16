#include "evpfft.h"
#include "defgrad_dcmp.h"
#include "euler.h"
#include "XdmfUniformGridWriter.h"
#include "math_functions.h"
#include <stdio.h>
#include <fstream>
#include <sys/stat.h>
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

void calc_IPF_colors(
  const MatrixTypeRealDual& ag,
  const MatrixTypeRealDual& defgrade,
  DFMatrixKokkos <unsigned int> &elem_IPF_colors)
{
  // This version of calc_IPF_colors is specific to LS-EVPFFT
  // because the orientation is updated inside here

  size_t npts1 = ag.dims(3);
  size_t npts2 = ag.dims(4);
  size_t npts3 = ag.dims(5);
  MatrixTypeRealDevice eulers(3);

  FOR_ALL(
        k, 1, npts3+1,
        j, 1, npts2+1,
        i, 1, npts1+1, {
      double pi = 3.141592653589793;
      real_t aa_[3*3];
      real_t Fe_local_[3*3];
      real_t V_[3*3]; 
      real_t R_[3*3];
      ViewMatrixTypeReal aa(aa_,3,3);
      ViewMatrixTypeReal Fe_local(Fe_local_,3,3);
      ViewMatrixTypeReal V(V_,3,3);
      ViewMatrixTypeReal R(R_,3,3);

      for (int ii = 1; ii <= 3; ii++) {
        for (int jj = 1; jj <= 3; jj++) {
          Fe_local(ii,jj) = defgrade(ii,jj,i,j,k);
        }
      }

      defgrad_dcmp(Fe_local.pointer(), V.pointer(), R.pointer());

      for (int ii = 1; ii <= 3; ii++) {
        for (int jj = 1; jj <= 3; jj++) {
          aa(ii,jj) = 0.0;
          for (int kk = 1; kk <= 3; kk++) {
            aa(ii,jj) =+ ag(ii,kk,i,j,k) * R(kk,jj);
          }
        }
      }
      euler(1, eulers(1), eulers(2), eulers(3), aa.pointer());
  
      double hkl[3]; //set reference direction for IPF
      hkl[0] = 0.;
      hkl[1] = 0.;
      hkl[2] = 1.;
      for (int ii=0;ii<3;ii++){
          hkl[ii] = hkl[ii]/sqrt(pow(hkl[0],2.)+pow(hkl[1],2.)+pow(hkl[2],2.));
      }

      //quatB convert Bunge Euler angles (radian) to quaternion
      double q[4];
      q[0] = sin(0.5*eulers(2))*cos(0.5*(eulers(1)-eulers(3)));
      q[1] = sin(0.5*eulers(2))*sin(0.5*(eulers(1)-eulers(3)));
      q[2] = cos(0.5*eulers(2))*sin(0.5*(eulers(1)+eulers(3)));
      q[3] = cos(0.5*eulers(2))*cos(0.5*(eulers(1)+eulers(3)));

      double dout[3];
      //invqrvec use quaterion to rotate vector in inverse sense, so is in reference to crystal system
      double t1 = pow(q[3],2.) - pow(q[0],2.) - pow(q[1],2.) - pow(q[2],2.);

      for (int ii=0;ii<3;ii++){
          dout[ii] = t1*hkl[ii];

          double tmp = 0.;
          for (int jj=0;jj<3;jj++){
              tmp = tmp + (q[jj]*hkl[jj]);
          }

          dout[ii] = dout[ii] + 2.*q[ii]*tmp;
          int kk;
          for (int jj=0;jj<3;jj++){
              if (ii!=jj){
                  for (int ijk=0;ijk<3;ijk++){
                      if (ii!=ijk and jj!=ijk){kk = ijk;}
                  }
                  int kl = jj-ii;
                  if (kl==2) {kl = -1;}
                  if (kl==-2) {kl = 1;}

                  dout[ii] = dout[ii] + 2.*q[kk]*q[3]*float(kl)*hkl[jj];
              }
          }
      }

      // Calculate RGB values from dout according to Groever's alg to get TSL color scheme
      dout[0] = fabs(dout[0]);
      dout[1] = fabs(dout[1]);
      dout[2] = fabs(dout[2]);
      double dmax = 0.;
      double dmin = 9999.;
      int bottom = 0;
      int middle = 0;
      int top = 0;

      for (int jj=0;jj<3;jj++){
          if (dmax < dout[jj]){
              dmax = dout[jj];
              top = jj;
          }
          if (dmin > dout[jj]){
              dmin = dout[jj];
              bottom = jj;
          }
      }
      for (int jj=0;jj<3;jj++){
          if (jj != bottom and jj != top){middle = jj;}
      }

      if (dout[top] > 1.){dout[top] = 1.;}
      double phi = fmin(pi/4.,acos(dout[top]));
      double thett = atan2(dout[bottom],dout[middle]);

      double red;
      double blue;
      double green;
      red = 1. - sin(1.6 * phi);
      green = (1. - red) * cos(2.*thett);
      blue = (1. - red) * sin(2.*thett);

      double drmax = fmax(red,fmax(green,blue));
      double correct = 1./drmax;
      red = red*correct;
      green = green*correct;
      blue = blue*correct;
      
      elem_IPF_colors(1,i,j,k) = (unsigned int) 255.*red;
      elem_IPF_colors(2,i,j,k) = (unsigned int) 255.*green;
      elem_IPF_colors(3,i,j,k) = (unsigned int) 255.*blue;

  });
  Kokkos::fence();
  elem_IPF_colors.update_host();
} 

void EVPFFT::write_micro_state_xdmf()
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

  // write jgrain
  write_N1xN2xN3 <int> (jgrain.pointer(), npts1, npts2, npts3, "grain_id",
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
void EVPFFT::write_micro_state_pvtu()
{
  Profiler profiler(__FUNCTION__);

  std::string tmp_str;
  std::stringstream mixed_str;

  char name[128];
  char pfilename[128];
  char filename[128];
  char pvtu_dir[128];
  char dirname[128];
  char subdirname[128];
  char tmp[128];

  int num_dims = 3; //EVPFFT always 3
  int graphics_idx = imicro;
  int num_nodes_in_elem = 8;
  int num_points = (npts1+1)*(npts2+1)*(npts3+1);
  int num_cells = npts1*npts2*npts3;

  MPI_Barrier(MPI_COMM_WORLD);

  const size_t n = 3*npts1_g*npts2_g*npts3_g; // for detFavg (1), defgradinvavgc (9), and velgradavg (9)

  // calculate elastic-strain
  MatrixTypeRealDual eel (3,3,npts1,npts2,npts3);
  calculate_eel(eel);

  // Some arrays for results
  MatrixTypeRealDual sym_disgrad (3,3,npts1,npts2,npts3);
  MatrixTypeRealDual sym_edotp (3,3,npts1,npts2,npts3);
  MatrixTypeRealDual sym_eel (3,3,npts1,npts2,npts3);
  DFMatrixKokkos <unsigned int> IPF_colors (3,npts1,npts2,npts3);

  // Symmetrize some fields
  symmetrize(disgrad, sym_disgrad);
  symmetrize(edotp, sym_edotp);
  symmetrize(eel, sym_eel);

  // Calculate IPF colors relative to z axis
  calc_IPF_colors(ag,defgrade,IPF_colors);


  sg.update_host();
  edotp.update_host();

  // Calculate point positions
  MatrixTypeRealDual defgradavg_dual(3,3);
  MatrixTypeRealDual uf(3,npts1_g,npts2_g,npts3_g);
  MatrixTypeRealDual ufintp(3,npts1+1,npts2+1,npts3+1);
  MatrixTypeRealDual xintp(3,npts1+1,npts2+1,npts3+1);
  MatrixTypeIntHost pid(npts1+1,npts2+1,npts3+1);

  for (int ii = 1; ii <= 3; ii++) {
    for (int jj = 1; jj <= 3; jj++) {
      defgradavg_dual.host(ii,jj) = defgradavg(ii,jj);
    }
  }
  defgradavg_dual.update_device();
  
  FOR_ALL_CLASS(
          kz, 1, npts3+1,
          ky, 1, npts2+1,
          kx, 1, npts1+1, {

        double xtmp[3];
        xtmp[0] = double(kx);
        xtmp[1] = double(ky);
        xtmp[2] = double(kz);
        for (int ii = 1; ii <= 3; ii++) {
          real_t dum = 0.0;
          for (int jj = 1; jj <= 3; jj++) {
            dum += defgradavg_dual(ii,jj)*xtmp[jj-1];
          }
          uf(ii,kx+local_start1,ky+local_start2,kz+local_start3) = xnode(ii,kx,ky,kz) - dum;
        }
  }); // end FOR_ALL_CLASS
  Kokkos::fence();
  uf.update_host();

  MPI_Allreduce(MPI_IN_PLACE, uf.host_pointer(), n, MPI_REAL_T, MPI_SUM, mpi_comm);

  FOR_ALL_CLASS(
          kz, 1, npts3+2,
          ky, 1, npts2+2,
          kx, 1, npts1+2, {

        int kxn1;
        int kxn2;
        int kyn1;
        int kyn2;
        int kzn1;
        int kzn2;
        // wrapped ids of neighbring points
        if ((kx+local_start1)==1) {
          kxn1 = kx+local_start1;
          kxn2 = npts1_g;
        } else if ((kx+local_start1)==(npts1_g+1)) {
          kxn1 = 1;
          kxn2 = npts1_g;
        } else {
          kxn1 = kx+local_start1;
          kxn2 = kx+local_start1-1;
        }
        if ((ky+local_start2)==1) {
          kyn1 = ky+local_start2;
          kyn2 = npts2_g;
        } else if ((ky+local_start2)==(npts2_g+1)) {
          kyn1 = 1;
          kyn2 = npts2_g;
        } else {
          kyn1 = ky+local_start2;
          kyn2 = ky+local_start2-1;
        }
        if ((kz+local_start3)==1) {
          kzn1 = kz+local_start3;
          kzn2 = npts3_g;
        } else if ((kz+local_start3)==(npts3_g+1)) {
          kzn1 = 1;
          kzn2 = npts3_g;
        } else {
          kzn1 = kz+local_start3;
          kzn2 = kz+local_start3-1;
        }

        for (int ii = 1; ii <= 3; ii++) {
          ufintp(ii,kx,ky,kz) = 0.125*(uf(ii,kxn1,kyn1,kzn1)+uf(ii,kxn2,kyn1,kzn1)+uf(ii,kxn1,kyn2,kzn1)+uf(ii,kxn2,kyn2,kzn1)+
          uf(ii,kxn1,kyn1,kzn2)+uf(ii,kxn2,kyn1,kzn2)+uf(ii,kxn1,kyn2,kzn2)+uf(ii,kxn2,kyn2,kzn2) );
        }
  }); // end FOR_ALL_CLASS
  Kokkos::fence();
  ufintp.update_host();

  FOR_ALL_CLASS(
          kz, 1, npts3+2,
          ky, 1, npts2+2,
          kx, 1, npts1+2, {

        double xtmp[3];
        xtmp[0] = double(kx+local_start1)-0.5;
        xtmp[1] = double(ky+local_start2)-0.5;
        xtmp[2] = double(kz+local_start3)-0.5;

        for (int ii = 1; ii <= 3; ii++) {
          real_t dum = 0.0;
          for (int jj = 1; jj <= 3; jj++) {
            dum += defgradavg_dual(ii,jj)*xtmp[jj-1];
          }
          xintp(ii,kx,ky,kz) = (dum + ufintp(ii,kx,ky,kz))*delt(ii);
        }
  }); // end FOR_ALL_CLASS
  Kokkos::fence();
  xintp.update_host();

  sprintf(name,"MicroState");
  sprintf(pvtu_dir, "pvtu");
  sprintf(subdirname, "%s_%05lu", name, graphics_idx);

  // mkdir if needed
  struct stat st;
  if (my_rank == 0) {
    if (stat(pvtu_dir, &st) != 0) {
      // str_stream.str("");
      // str_stream << "mkdir" << " " << pvtu_dir;
      // system(str_stream.str().c_str());
      //tmp = "";
      sprintf(tmp, "mkdir %s",pvtu_dir);
      system(tmp);
    }
    sprintf(dirname,"%s/%s",pvtu_dir,subdirname);
    if (stat(dirname, &st) != 0) {
      // str_stream.str("");
      // str_stream << "mkdir" << " " << vtu_dir;
      // system(str_stream.str().c_str());
      //tmp = "";
      sprintf(tmp, "mkdir %s", dirname);
      system(tmp);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  //  ---------------------------------------------------------------------------
  //  Write the PVTU file (only done by a single rank)
  //  ---------------------------------------------------------------------------
  unsigned long int byte_offset = 0;
  int block_header_size = sizeof(unsigned long int);
  unsigned long int data_block_size;

  if (my_rank == 0){
    sprintf(pfilename, "%s/%s_%05lu.pvtu",pvtu_dir, name, graphics_idx);
    std::ofstream pout;  // FILE *out;
    pout.open(pfilename,std::ofstream::binary);

    //  Write Header
    tmp_str = "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    pout.write(tmp_str.c_str(), tmp_str.length());
    tmp_str = "<PUnstructuredGrid GhostLevel=\"1\">\n"; //unsure of exact correct usage of ghostlevel, leaving as 1 for now
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  Write Field Data (only part in pvtu using byte_offset)
    tmp_str = "<FieldData>\n";
    pout.write(tmp_str.c_str(), tmp_str.length());
    mixed_str.str("");
    mixed_str << "<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
    byte_offset += sizeof(double)+block_header_size;
    tmp_str = mixed_str.str();
    pout.write(tmp_str.c_str(), tmp_str.length());
    tmp_str = "</FieldData>\n";
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  Write PPoints Section
    tmp_str = "<PPoints>\n";
    pout.write(tmp_str.c_str(), tmp_str.length());
    mixed_str.str("");
    mixed_str << "<PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>\n";
    tmp_str = mixed_str.str();
    pout.write(tmp_str.c_str(), tmp_str.length());
    tmp_str = "</PPoints>\n";
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  Write PPoint Data List
    tmp_str = "<PPointData>\n";
    pout.write(tmp_str.c_str(), tmp_str.length());
    //
    // //  Velocity
    // mixed_str.str("");
    // mixed_str << "<PDataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\"/>\n";
    // tmp_str = mixed_str.str();
    // pout.write(tmp_str.c_str(), tmp_str.length());
    //
    tmp_str = "</PPointData>\n";
    pout.write(tmp_str.c_str(), tmp_str.length());

    //  Write PCell Data List
    tmp_str = "<PCellData>\n";
    pout.write(tmp_str.c_str(), tmp_str.length());
    //
    //SCALARS
    //  Rank
    mixed_str.str("");
    mixed_str << "<PDataArray type=\"Int32\" Name=\"rank\"/>\n";
    tmp_str = mixed_str.str();
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  Phase ID
    mixed_str.str("");
    mixed_str << "<PDataArray type=\"Int32\" Name=\"phase_id\"/>\n";
    tmp_str = mixed_str.str();
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  Grain ID
    mixed_str.str("");
    mixed_str << "<PDataArray type=\"Int32\" Name=\"grain_id\"/>\n";
    tmp_str = mixed_str.str();
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  VM strain
    mixed_str.str("");
    mixed_str << "<PDataArray type=\"Float64\" Name=\"vm-strain\"/>\n";
    tmp_str = mixed_str.str();
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  VM stress
    mixed_str.str("");
    mixed_str << "<PDataArray type=\"Float64\" Name=\"vm-stress\"/>\n";
    tmp_str = mixed_str.str();
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  VM plastic strain rate
    mixed_str.str("");
    mixed_str << "<PDataArray type=\"Float64\" Name=\"vm-pl-strain-rate\"/>\n";
    tmp_str = mixed_str.str();
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  VM elastic strain
    mixed_str.str("");
    mixed_str << "<PDataArray type=\"Float64\" Name=\"vm-el-strain\"/>\n";
    tmp_str = mixed_str.str();
    pout.write(tmp_str.c_str(), tmp_str.length());
    //
    //VECTORS
    //  Strain
    mixed_str.str("");
    mixed_str << "<PDataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"6\"/>\n";
    tmp_str = mixed_str.str();
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  Stress
    mixed_str.str("");
    mixed_str << "<PDataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"6\"/>\n";
    tmp_str = mixed_str.str();
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  Plastic Strain Rate
    mixed_str.str("");
    mixed_str << "<PDataArray type=\"Float64\" Name=\"pl-strain-rate\" NumberOfComponents=\"6\"/>\n";
    tmp_str = mixed_str.str();
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  Elastic Strain
    mixed_str.str("");
    mixed_str << "<PDataArray type=\"Float64\" Name=\"el-strain\" NumberOfComponents=\"6\"/>\n";
    tmp_str = mixed_str.str();
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  IPFColor
    mixed_str.str("");
    mixed_str << "<PDataArray type=\"UInt8\" Name=\"IPFColor\" NumberOfComponents=\"3\"/>\n";
    tmp_str = mixed_str.str();
    pout.write(tmp_str.c_str(), tmp_str.length());

    // //  Stress
    // mixed_str.str("");
    // mixed_str << "<PDataArray type=\"Float64\" Name=\"Stress\" NumberOfComponents=\"9\"/>\n";
    // tmp_str = mixed_str.str();
    // pout.write(tmp_str.c_str(), tmp_str.length());
    // //  Density
    // mixed_str.str("");
    // mixed_str << "<PDataArray type=\"Float64\" Name=\"density\"/>\n";
    // tmp_str = mixed_str.str();
    // pout.write(tmp_str.c_str(), tmp_str.length());

    //
    tmp_str = "</PCellData>\n";
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  Write Partition Piece List
    char piecename[128];
    for (int rank = 0; rank < num_ranks; rank++){
      sprintf(piecename, "%s/%s_%05lu_%05lu.vtu",subdirname, name, graphics_idx, rank);
      mixed_str.str("");
      mixed_str << "<Piece Source=\"" << piecename << "\"/>\n";
      tmp_str = mixed_str.str();
      pout.write(tmp_str.c_str(), tmp_str.length());
    }
    //
    tmp_str = "</PUnstructuredGrid>\n";
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  Write Appended Data
    tmp_str = "<AppendedData encoding=\"raw\">\n_";
    pout.write(tmp_str.c_str(), tmp_str.length());
    //  Write Time Value
    data_block_size = sizeof(double);
    pout.write((char *) &data_block_size,block_header_size);
    double  time_value = double(imicro)*tdot;
    pout.write((char *) &time_value,sizeof(time_value));
    pout.put('\n');
    tmp_str = "</AppendedData>\n";
    pout.write(tmp_str.c_str(), tmp_str.length());  

    tmp_str = "</VTKFile>\n";
    pout.write(tmp_str.c_str(), tmp_str.length());

    pout.close();
  }

  //  ---------------------------------------------------------------------------
  //  Write the VTU file
  //  ---------------------------------------------------------------------------
  sprintf(filename, "%s/%s/%s_%05lu_%05lu.vtu",pvtu_dir ,subdirname, name, graphics_idx, my_rank);
  // filename has the full string
  
  std::ofstream out;  // FILE *out;
  out.open(filename,std::ofstream::binary);

  byte_offset = 0;
  
  //  Write Header
  //tmp_str = "<?xml version=\"1.0\"?>\n";
  //out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "<UnstructuredGrid>\n";
  out.write(tmp_str.c_str(), tmp_str.length());

  //  Write Time Value Header
  tmp_str = "<FieldData>\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  mixed_str.str("");
  mixed_str << "<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
  byte_offset += sizeof(double)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "</FieldData>\n";
  out.write(tmp_str.c_str(), tmp_str.length());

  mixed_str.str("");
  mixed_str << "<Piece NumberOfPoints=\"" << num_points << "\" NumberOfCells=\"" << num_cells << "\">\n";
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());

  //  **Write Points Header**
  tmp_str = "<Points>\n";
  out.write(tmp_str.c_str(), tmp_str.length());   
  //tmp_str = "<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\">\n"; 
  mixed_str.str("");
  mixed_str << "<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
  byte_offset += num_points*num_dims*sizeof(double)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "</DataArray>\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  // tmp_str = "</DataArray>\n";
  // out.write(tmp_str.c_str(), tmp_str.length());  
  tmp_str = "</Points>\n";
  out.write(tmp_str.c_str(), tmp_str.length());

  //  **Write Cells Header**
  tmp_str = "<Cells>\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  //  Connectivity
  mixed_str.str("");
  mixed_str << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
  byte_offset += num_cells*num_nodes_in_elem*sizeof(int)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  //  Offsets
  mixed_str.str("");
  mixed_str << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n";
  byte_offset += num_cells*sizeof(int)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  //  Types
  mixed_str.str("");
  mixed_str << "<DataArray type=\"Int32\" Name=\"types\" format=\"appended\" offset=\"" << byte_offset << "\"/>\n"; 
  byte_offset += num_cells*sizeof(int)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "</Cells>\n";
  out.write(tmp_str.c_str(), tmp_str.length());

  //  **Write Point Data Headers**
  tmp_str = "<PointData>\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  //
  // //  Velocity
  // mixed_str.str("");
  // mixed_str << "<DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
  // byte_offset += num_points*num_dims*sizeof(double)+block_header_size;
  // tmp_str = mixed_str.str();
  // out.write(tmp_str.c_str(), tmp_str.length());
  // tmp_str = "</DataArray>\n";
  // out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "</PointData>\n";
  out.write(tmp_str.c_str(), tmp_str.length());

  //  **Write Cell Data Headers**
  tmp_str = "<CellData>\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  //
  //SCALARS
  //  Rank
  mixed_str.str("");
  mixed_str << "<DataArray type=\"Int32\" Name=\"rank\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
  byte_offset += num_cells*sizeof(int)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "</DataArray>\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  //  Phase ID
  mixed_str.str("");
  mixed_str << "<DataArray type=\"Int32\" Name=\"phase_id\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
  byte_offset += num_cells*sizeof(int)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "</DataArray>\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  //  Grain ID
  mixed_str.str("");
  mixed_str << "<DataArray type=\"Int32\" Name=\"grain_id\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
  byte_offset += num_cells*sizeof(int)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "</DataArray>\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  //  VM Strain
  mixed_str.str("");
  mixed_str << "<DataArray type=\"Float64\" Name=\"vm-strain\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
  byte_offset += num_cells*sizeof(double)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "</DataArray>\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  //  VM Stress
  mixed_str.str("");
  mixed_str << "<DataArray type=\"Float64\" Name=\"vm-stress\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
  byte_offset += num_cells*sizeof(double)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "</DataArray>\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  //  VM Plastic Strain Rate
  mixed_str.str("");
  mixed_str << "<DataArray type=\"Float64\" Name=\"vm-pl-strain-rate\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
  byte_offset += num_cells*sizeof(double)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "</DataArray>\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  //  VM Elastic Strain
  mixed_str.str("");
  mixed_str << "<DataArray type=\"Float64\" Name=\"vm-el-strain\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
  byte_offset += num_cells*sizeof(double)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "</DataArray>\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  //
  //VECTORS
  //  Strain
  mixed_str.str("");
  mixed_str << "<DataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"6\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
  byte_offset += 6*num_cells*sizeof(double)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "</DataArray>\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  //  Stress
  mixed_str.str("");
  mixed_str << "<DataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"6\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
  byte_offset += 6*num_cells*sizeof(double)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "</DataArray>\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  //  Plastic Strain Rate
  mixed_str.str("");
  mixed_str << "<DataArray type=\"Float64\" Name=\"pl-strain-rate\" NumberOfComponents=\"6\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
  byte_offset += 6*num_cells*sizeof(double)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "</DataArray>\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  //  Elastic Strain
  mixed_str.str("");
  mixed_str << "<DataArray type=\"Float64\" Name=\"el-strain\" NumberOfComponents=\"6\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
  byte_offset += 6*num_cells*sizeof(double)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "</DataArray>\n";
  out.write(tmp_str.c_str(), tmp_str.length());
  //  IPFColor
  mixed_str.str("");
  mixed_str << "<DataArray type=\"UInt8\" Name=\"IPFColor\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
  byte_offset += num_cells*num_dims*sizeof(unsigned char)+block_header_size;
  tmp_str = mixed_str.str();
  out.write(tmp_str.c_str(), tmp_str.length());
  tmp_str = "</DataArray>\n";
  out.write(tmp_str.c_str(), tmp_str.length());

  // //  Stress
  // mixed_str.str("");
  // mixed_str << "<DataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"9\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
  // byte_offset += num_cells*sizeof(double)+block_header_size;
  // tmp_str = mixed_str.str();
  // out.write(tmp_str.c_str(), tmp_str.length());
  // tmp_str = "</DataArray>\n";
  // out.write(tmp_str.c_str(), tmp_str.length());
  // //  Density
  // mixed_str.str("");
  // mixed_str << "<DataArray type=\"Float64\" Name=\"density\" format=\"appended\" offset=\"" << byte_offset << "\">\n";
  // byte_offset += num_cells*sizeof(double)+block_header_size;
  // tmp_str = mixed_str.str();
  // out.write(tmp_str.c_str(), tmp_str.length());
  // tmp_str = "</DataArray>\n";
  // out.write(tmp_str.c_str(), tmp_str.length());
  //
  tmp_str = "</CellData>\n";
  out.write(tmp_str.c_str(), tmp_str.length());

  //  Write Mesh Close
  tmp_str = "</Piece>\n";
  out.write(tmp_str.c_str(), tmp_str.length());  
  tmp_str = "</UnstructuredGrid>\n";
  out.write(tmp_str.c_str(), tmp_str.length());  

  //  Write Appended Data
  tmp_str = "<AppendedData encoding=\"raw\">\n_";
  out.write(tmp_str.c_str(), tmp_str.length());  

  //  Write Time Value
  data_block_size = sizeof(double);
  out.write((char *) &data_block_size,block_header_size);
  double  time_value = double(imicro)*tdot;
  out.write((char *) &time_value,sizeof(time_value));

  //  **Write Points Data**
  //double coords[num_points*num_dims];
  data_block_size = num_points*num_dims*sizeof(double);
  out.write((char *) &data_block_size,block_header_size);
  double coord_tmp;
  // for (int node_gid = 0; node_gid < num_points; node_gid++){
  //   for (int dim = 0; dim < num_dims; dim++){
  //     coord_tmp = node_coords(node_gid,dim);
  //     out.write((char *) &coord_tmp,sizeof(coord_tmp));
  //   }
  // }
  // Change to looping over grid indicies
  int ic = -1;
  for (int kz = 1; kz <= npts3+1; kz++){
    for (int ky = 1; ky <= npts2+1; ky++){
      for (int kx = 1; kx <= npts1+1; kx++){
        ic += 1;
        for (int dim = 1; dim <= num_dims; dim++){
          coord_tmp = xintp.host(dim,kx,ky,kz);
          out.write((char *) &coord_tmp,sizeof(coord_tmp));
        }
        pid(kx,ky,kz) = ic;
      }
    }
  }

  //  **Write Cells Data**
  //int connect[num_cells*mesh.num_nodes_in_elem];
  data_block_size = num_cells*num_nodes_in_elem*sizeof(int);
  out.write((char *) &data_block_size,block_header_size);
  int connect_tmp;
  // for (int elem_gid = 0; elem_gid < num_cells; elem_gid++){
  //   for (int node_lid = 0; node_lid < num_nodes_in_elem; node_lid++){
  //     connect_tmp = nodes_in_elem(elem_gid, node_lid);
  //     out.write((char *) &connect_tmp,sizeof(connect_tmp));
  //   }
  // }
  // Change to looping over grid indicies
  for (int kz = 1; kz <= npts3; kz++){
    for (int ky = 1; ky <= npts2; ky++){
      for (int kx = 1; kx <= npts1; kx++){
        // Cumbersome grid-based listing of cell points
        connect_tmp = pid(kx,ky,kz);
        out.write((char *) &connect_tmp,sizeof(connect_tmp));
        connect_tmp = pid(kx+1,ky,kz);
        out.write((char *) &connect_tmp,sizeof(connect_tmp));
        connect_tmp = pid(kx+1,ky+1,kz);
        out.write((char *) &connect_tmp,sizeof(connect_tmp));
        connect_tmp = pid(kx,ky+1,kz);
        out.write((char *) &connect_tmp,sizeof(connect_tmp));
        connect_tmp = pid(kx,ky,kz+1);
        out.write((char *) &connect_tmp,sizeof(connect_tmp));
        connect_tmp = pid(kx+1,ky,kz+1);
        out.write((char *) &connect_tmp,sizeof(connect_tmp));
        connect_tmp = pid(kx+1,ky+1,kz+1);
        out.write((char *) &connect_tmp,sizeof(connect_tmp));
        connect_tmp = pid(kx,ky+1,kz+1);
        out.write((char *) &connect_tmp,sizeof(connect_tmp));
      }
    }
  }
  //int offsets[num_cells];
  data_block_size = num_cells*sizeof(int);
  out.write((char *) &data_block_size,block_header_size);
  int offset_tmp;
  for (int elem_gid = 0; elem_gid < num_cells; elem_gid++){
    offset_tmp = (elem_gid+1)*num_nodes_in_elem;
    out.write((char *) &offset_tmp,sizeof(offset_tmp));
  }
  //int types[num_cells];
  data_block_size = num_cells*sizeof(int);
  out.write((char *) &data_block_size,block_header_size);
  int type_tmp;
  for (int elem_gid = 0; elem_gid < num_cells; elem_gid++){
    type_tmp = 12;
    out.write((char *) &type_tmp,sizeof(type_tmp));
  }

  //  **Write Point Data**
  //SCALARS float
  // //  Velocity
  // data_block_size = num_points*num_dims*sizeof(double);
  // out.write((char *) &data_block_size,block_header_size);
  // for (int node_gid = 0; node_gid < num_points; node_gid++){
  //   for (int dim = 0; dim < num_dims; dim++){
  //     //out.write((char *) &node_vel.host(1, node_gid, dim),sizeof(double));
  //   }
  // }

  //  **Write Cell Data**
  //SCALARS
  //  Rank
  data_block_size = num_cells*sizeof(int);
  out.write((char *) &data_block_size,block_header_size);
  for (int elem_gid = 0; elem_gid < num_cells; elem_gid++){
      out.write((char *) &my_rank,sizeof(int));
  }
  //  Phase ID
  data_block_size = num_cells*sizeof(int);
  out.write((char *) &data_block_size,block_header_size);
  for (int kz = 1; kz <= npts3; kz++){
    for (int ky = 1; ky <= npts2; ky++){
      for (int kx = 1; kx <= npts1; kx++){
        out.write((char *) &jphase.host(kx,ky,kz),sizeof(int));
      }
    }
  }
  //  Grain ID
  data_block_size = num_cells*sizeof(int);
  out.write((char *) &data_block_size,block_header_size);
  for (int kz = 1; kz <= npts3; kz++){
    for (int ky = 1; ky <= npts2; ky++){
      for (int kx = 1; kx <= npts1; kx++){
        out.write((char *) &jgrain(kx,ky,kz),sizeof(int));
      }
    }
  }
  // for writing vm fields 
  MatrixTypeRealDual vm_field (npts1,npts2,npts3);
  //  VM strain
  calculate_vm_field(sym_disgrad, vm_field, 1);
  data_block_size = num_cells*sizeof(double);
  out.write((char *) &data_block_size,block_header_size);
  for (int kz = 1; kz <= npts3; kz++){
    for (int ky = 1; ky <= npts2; ky++){
      for (int kx = 1; kx <= npts1; kx++){
        out.write((char *) &vm_field.host(kx,ky,kz),sizeof(double));
      }
    }
  }
  //  VM stress
  calculate_vm_field(sg, vm_field, 1);
  data_block_size = num_cells*sizeof(double);
  out.write((char *) &data_block_size,block_header_size);
  for (int kz = 1; kz <= npts3; kz++){
    for (int ky = 1; ky <= npts2; ky++){
      for (int kx = 1; kx <= npts1; kx++){
        out.write((char *) &vm_field.host(kx,ky,kz),sizeof(double));
      }
    }
  }
  //  VM plastic strain rate
  calculate_vm_field(sym_edotp, vm_field, 1);
  data_block_size = num_cells*sizeof(double);
  out.write((char *) &data_block_size,block_header_size);
  for (int kz = 1; kz <= npts3; kz++){
    for (int ky = 1; ky <= npts2; ky++){
      for (int kx = 1; kx <= npts1; kx++){
        out.write((char *) &vm_field.host(kx,ky,kz),sizeof(double));
      }
    }
  }
  //  VM elastic strain
  calculate_vm_field(sym_eel, vm_field, 1);
  data_block_size = num_cells*sizeof(double);
  out.write((char *) &data_block_size,block_header_size);
  for (int kz = 1; kz <= npts3; kz++){
    for (int ky = 1; ky <= npts2; ky++){
      for (int kx = 1; kx <= npts1; kx++){
        out.write((char *) &vm_field.host(kx,ky,kz),sizeof(double));
      }
    }
  }
  //
  //Vectors
  //  Strain
  data_block_size = 6*num_cells*sizeof(double);
  out.write((char *) &data_block_size,block_header_size);
  for (int kz = 1; kz <= npts3; kz++){
    for (int ky = 1; ky <= npts2; ky++){
      for (int kx = 1; kx <= npts1; kx++){
        // Manually write out order until a better way is found
        out.write((char *) &sym_disgrad.host(1,1,kx,ky,kz),sizeof(double));
        out.write((char *) &sym_disgrad.host(2,2,kx,ky,kz),sizeof(double));
        out.write((char *) &sym_disgrad.host(3,3,kx,ky,kz),sizeof(double));
        out.write((char *) &sym_disgrad.host(1,2,kx,ky,kz),sizeof(double));
        out.write((char *) &sym_disgrad.host(2,3,kx,ky,kz),sizeof(double));
        out.write((char *) &sym_disgrad.host(1,3,kx,ky,kz),sizeof(double));
      }
    }
  }
  //  Stress
  data_block_size = 6*num_cells*sizeof(double);
  out.write((char *) &data_block_size,block_header_size);
  for (int kz = 1; kz <= npts3; kz++){
    for (int ky = 1; ky <= npts2; ky++){
      for (int kx = 1; kx <= npts1; kx++){
        // Manually write out order until a better way is found
        out.write((char *) &sg.host(1,1,kx,ky,kz),sizeof(double));
        out.write((char *) &sg.host(2,2,kx,ky,kz),sizeof(double));
        out.write((char *) &sg.host(3,3,kx,ky,kz),sizeof(double));
        out.write((char *) &sg.host(1,2,kx,ky,kz),sizeof(double));
        out.write((char *) &sg.host(2,3,kx,ky,kz),sizeof(double));
        out.write((char *) &sg.host(1,3,kx,ky,kz),sizeof(double));
      }
    }
  }
  //  Plastic Strain rate
  data_block_size = 6*num_cells*sizeof(double);
  out.write((char *) &data_block_size,block_header_size);
  for (int kz = 1; kz <= npts3; kz++){
    for (int ky = 1; ky <= npts2; ky++){
      for (int kx = 1; kx <= npts1; kx++){
        // Manually write out order until a better way is found
        out.write((char *) &sym_edotp.host(1,1,kx,ky,kz),sizeof(double));
        out.write((char *) &sym_edotp.host(2,2,kx,ky,kz),sizeof(double));
        out.write((char *) &sym_edotp.host(3,3,kx,ky,kz),sizeof(double));
        out.write((char *) &sym_edotp.host(1,2,kx,ky,kz),sizeof(double));
        out.write((char *) &sym_edotp.host(2,3,kx,ky,kz),sizeof(double));
        out.write((char *) &sym_edotp.host(1,3,kx,ky,kz),sizeof(double));
      }
    }
  }
  //  Elastic Strain
  data_block_size = 6*num_cells*sizeof(double);
  out.write((char *) &data_block_size,block_header_size);
  for (int kz = 1; kz <= npts3; kz++){
    for (int ky = 1; ky <= npts2; ky++){
      for (int kx = 1; kx <= npts1; kx++){
        // Manually write out order until a better way is found
        out.write((char *) &sym_eel.host(1,1,kx,ky,kz),sizeof(double));
        out.write((char *) &sym_eel.host(2,2,kx,ky,kz),sizeof(double));
        out.write((char *) &sym_eel.host(3,3,kx,ky,kz),sizeof(double));
        out.write((char *) &sym_eel.host(1,2,kx,ky,kz),sizeof(double));
        out.write((char *) &sym_eel.host(2,3,kx,ky,kz),sizeof(double));
        out.write((char *) &sym_eel.host(1,3,kx,ky,kz),sizeof(double));
      }
    }
  }
  //  IPFColor
  data_block_size = num_cells*num_dims*sizeof(unsigned char);
  out.write((char *) &data_block_size,block_header_size);
  for (int kz = 1; kz <= npts3; kz++){
    for (int ky = 1; ky <= npts2; ky++){
      for (int kx = 1; kx <= npts1; kx++){
        for (int dim = 1; dim <= num_dims; dim++){
            out.write((char *) &IPF_colors.host(dim,kx,ky,kz),sizeof(unsigned char));
        }
      }
    }
  }

  // //  Stress
  // data_block_size = num_cells*sizeof(double);
  // out.write((char *) &data_block_size,block_header_size);
  // for (int elem_gid = 0; elem_gid < num_cells; elem_gid++){
  //   for (int ii = 0; ii < 9; ii++){
  //     //out.write((char *) &elem_stress.host(1,elem_gid,ii),sizeof(double));
  //   }
  // }
  // //  Density
  // data_block_size = num_cells*sizeof(double);
  // out.write((char *) &data_block_size,block_header_size);
  // for (int elem_gid = 0; elem_gid < num_cells; elem_gid++){
  //     out.write((char *) &elem_den.host(elem_gid),sizeof(double));
  // }

  
  out.put('\n');
  tmp_str = "</AppendedData>\n";
  out.write(tmp_str.c_str(), tmp_str.length());  

  //  Write File Close
  tmp_str = "</VTKFile>\n";
  out.write(tmp_str.c_str(), tmp_str.length());  

  out.close();

}

