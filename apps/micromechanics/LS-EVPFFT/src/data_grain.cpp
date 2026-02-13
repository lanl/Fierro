#include <stdio.h>
#include "evpfft.h"
#include "definitions.h"
#include "utilities.h"
#include "euler.h"
#include "math_functions.h"
#include "voigt.h"
#include "hdf5_io_functions.h"

using namespace utils;

int get_dataset_rank(hid_t hdf5_file_id, const char *dataset_name);

void EVPFFT::data_grain(const std::string & filetext)
{

  switch (cmd.micro_filetype)
  {
    case 0:
    read_classic_los_alamos_texture_file(filetext);
    break;

    case 1:
    read_hdf5_texture_file(filetext);
    break;

    case 2:
    read_vtk_lattice_structure(filetext);
    break;

    default:
      throw std::runtime_error("invalid micro_filetype.");
  }

#if 0
  // for debug
  mpi_io_real_t->write_ascii("debug_ph.txt", ph_array.host_pointer(), "");
  mpi_io_real_t->write_ascii("debug_th.txt", th_array.host_pointer(), "");
  mpi_io_real_t->write_ascii("debug_om.txt", om_array.host_pointer(), "");
  mpi_io_int->write_ascii("debug_jgr.txt", jgrain.host_pointer(), "");
  mpi_io_int->write_ascii("debug_jph.txt", jphase.host_pointer(), "");
#endif

  int nph1 = 0;
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

      if (jph == 1) nph1 = nph1 + 1;
      if (igas.host(jph) == 1) continue; // nothing needs to be done if gas phase
           
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
          c066.host(i,j) += caux66(i,j)*wgt;  
        }
      }

      } // end for ii
    } // end for jj
  } // end for kk

  // update device
  cg66.update_device();
  ag.update_device();
  Kokkos::fence();

  // Allreduce on c066
  MPI_Allreduce(MPI_IN_PLACE, c066.host_pointer(), c066.size(), MPI_REAL_T, MPI_SUM, mpi_comm);
  c066.update_device();
  
  // for single crystal only. it leads to fast convergence of the outer while loop for stress and strain fields
  // for (int i = 1; i <= 6; i++) {
  //   for (int j = 1; j <= 6; j++) {
  //     c066.host(i,j) *= 1000000.0;
  //   }
  // }
  // c066.update_device();

#if 0
  // debug print
  print_array_to_file(c066.host_pointer(), c066.size(), my_rank, "c066.txt");
#endif

  wph1 = (1.0*nph1)/(npts1_g*npts2_g*npts3_g);

  MatrixTypeRealHost s066(6,6);
  for (int i = 1; i <= 6; i++) {
    for (int j = 1; j <= 6; j++) {
      s066(i,j) = c066.host(i,j);
    }
  }
  invert_matrix <6> (s066.pointer());

  cb.chg_basis_3(c066.host_pointer(), c0.host_pointer(), 3, 6, cb.B_basis_host_pointer());
  cb.chg_basis_3(s066.pointer(), s0.host_pointer(), 3, 6, cb.B_basis_host_pointer());

  // update device
  c0.update_device();
  s0.update_device();
  Kokkos::fence();
}


void EVPFFT::read_classic_los_alamos_texture_file(const std::string & filetext)
{
  std::ifstream ur2;
  ur2.open(filetext);
  check_that_file_is_open(ur2, filetext.c_str());

  real_t ph, th, om;
  int iiii, jjjj, kkkk, jgr, jph;
  int i,j,k;

  for (int kk = 0; kk < npts3_g; kk++) {
    for (int jj = 0; jj < npts2_g; jj++) {
      for (int ii = 0; ii < npts1_g; ii++) {
      
        ur2 >> ph >> th >> om >> iiii >> jjjj >> kkkk >> jgr >> jph; CLEAR_LINE(ur2);
        if ( (ii >= local_start1 and ii <= local_end1) and
             (jj >= local_start2 and jj <= local_end2) and
             (kk >= local_start3 and kk <= local_end3) )
        {
          i = ii - local_start1 + 1;
          j = jj - local_start2 + 1;
          k = kk - local_start3 + 1;
          ph_array(i,j,k) = ph;
          th_array(i,j,k) = th;
          om_array(i,j,k) = om;
          jgrain(i,j,k)   = jgr;
          jphase.host(i,j,k)   = jph;
        }
        
      }
    }
  }

  // update device
  jphase.update_device();
  Kokkos::fence();

  ur2.close();
}

void EVPFFT::read_hdf5_texture_file(const std::string & filetext)
{
  // open hdf5 file
  hid_t hdf5_file_id = open_hdf5_file(filetext.c_str(), mpi_comm);
 
  // get all dataset ranks
  int EulerAnglesRank = get_dataset_rank(hdf5_file_id, cmd.EulerAnglesFullPath.c_str());
  int FeatureIdsRank = get_dataset_rank(hdf5_file_id, cmd.FeatureIdsFullPath.c_str());
  int PhaseIdsRank = get_dataset_rank(hdf5_file_id, cmd.PhaseIdsFullPath.c_str());

  if (EulerAnglesRank != 4) {
    throw std::runtime_error("EulerAngles should be stored as rank 4 array in hdf5 file");
  }
    
  // set dims and others
  // Note that hdf5 is stored in C layout, and when read from a c or c++ program
  // it reads in c-array-order but when read from fortran it reads in fortran-array-layout.
  // Therefore, the above arrays of dims should be reversed
  // Example: global_array_dims = [nx,ny,nz] should be [nz,ny,nx] 
  int global_array_dims[4] = {npts3_g,npts2_g,npts1_g,1};
  int EulerAngles_global_array_dims[4] = {npts3_g,npts2_g,npts1_g,3};
  int local_array_dims[4] = {npts3,npts2,npts1,1};
  int EulerAngles_local_array_dims[4] = {npts3,npts2,npts1,3};
  int start_coordinates[4] = {local_start3,local_start2,local_start1,0};
  int stride[4] = {1,1,1,1};
  
  // create filespace, memspace
  hid_t FeatureIds_filespace = create_hdf5_filespace(FeatureIdsRank, global_array_dims, local_array_dims,
                               start_coordinates, stride);
  hid_t PhaseIds_filespace = create_hdf5_filespace(PhaseIdsRank, global_array_dims, local_array_dims,
                             start_coordinates, stride);
  hid_t EulerAngles_filespace = create_hdf5_filespace(EulerAnglesRank, EulerAngles_global_array_dims,
                                EulerAngles_local_array_dims, start_coordinates, stride);
  hid_t FeatureIds_memspace = create_hdf5_memspace(FeatureIdsRank, local_array_dims, stride, 0);
  hid_t PhaseIds_memspace = create_hdf5_memspace(PhaseIdsRank, local_array_dims, stride, 0);
  hid_t EulerAngles_memspace = create_hdf5_memspace(EulerAnglesRank, EulerAngles_local_array_dims, stride, 0);
  
  // array to hold datasets
  FMatrix <real_t> EulerAngles (3,npts1,npts2,npts3);
  FMatrix <int> FeatureIds (1,npts1,npts2,npts3);
  FMatrix <int> PhaseIds (1,npts1,npts2,npts3);

  // read datasets
  read_hdf5_dataset <real_t> (cmd.EulerAnglesFullPath.c_str(), EulerAngles.pointer(), hdf5_file_id,
                              EulerAngles_memspace, EulerAngles_filespace, mpi_comm);
  read_hdf5_dataset <int> (cmd.FeatureIdsFullPath.c_str(), FeatureIds.pointer(), hdf5_file_id,
                           FeatureIds_memspace, FeatureIds_filespace, mpi_comm);
  read_hdf5_dataset <int> (cmd.PhaseIdsFullPath.c_str(), PhaseIds.pointer(), hdf5_file_id,
                           PhaseIds_memspace, PhaseIds_filespace, mpi_comm);
 
  // write data into the appropriate arrays
  real_t ph, th, om;
  int iiii, jjjj, kkkk, jgr, jph;
  int i,j,k;
  for (int kk = 1; kk <= npts3; kk++) {
    for (int jj = 1; jj <= npts2; jj++) {
      for (int ii = 1; ii <= npts1; ii++) {

        iiii = ii;
        jjjj = jj;
        kkkk = kk;
        ph = EulerAngles(1,ii,jj,kk);
        th = EulerAngles(2,ii,jj,kk);
        om = EulerAngles(3,ii,jj,kk);
        jgr = FeatureIds(1,ii,jj,kk); 
        jph = PhaseIds(1,ii,jj,kk);

        if (cmd.EulerAnglesUnit == 1) {
          // change from radian to degree
          ph = ph*180.0/pi;
          th = th*180.0/pi;
          om = om*180.0/pi;
        }

        ph_array(ii,jj,kk) = ph;
        th_array(ii,jj,kk) = th;
        om_array(ii,jj,kk) = om;
        jgrain(ii,jj,kk)   = jgr;
        jphase.host(ii,jj,kk) = jph;
        
      }
    }
  }

  // update device
  jphase.update_device();
  Kokkos::fence();
  
  H5Fclose(hdf5_file_id);
}

int get_dataset_rank(hid_t hdf5_file_id, const char *dataset_name)
{
  int dataset_rank;
  hid_t dataset = H5Dopen2(hdf5_file_id, dataset_name, H5P_DEFAULT);
  hid_t dataspace = H5Dget_space(dataset);
  dataset_rank = H5Sget_simple_extent_ndims(dataspace);
  H5Dclose(dataset);
  return dataset_rank;
}

void EVPFFT::read_vtk_lattice_structure(const std::string & filetext)
{
  std::ifstream ur2;
  ur2.open(filetext);
  check_that_file_is_open(ur2, filetext.c_str());

  // Skip header lines
  std::string line;
  while (std::getline(ur2, line)) {
    if (line == "LOOKUP_TABLE default") {
      break;
    }
  }

  int i,j,k;
  double grid_value;
  for (int kk = 0; kk < npts3_g; kk++) {
    for (int jj = 0; jj < npts2_g; jj++) {
      for (int ii = 0; ii < npts1_g; ii++) {
      
        ur2 >> grid_value;
        if ( (ii >= local_start1 and ii <= local_end1) and
             (jj >= local_start2 and jj <= local_end2) and
             (kk >= local_start3 and kk <= local_end3) )
        {
          i = ii - local_start1 + 1;
          j = jj - local_start2 + 1;
          k = kk - local_start3 + 1;
          ph_array(i,j,k)    = 0;
          th_array(i,j,k)    = 0;
          om_array(i,j,k)    = 0;
          jgrain(i,j,k)      = grid_value + 1;
          jphase.host(i,j,k) = grid_value + 1;
        }

      }
    }
  }
  // update device
  jphase.update_device();
  Kokkos::fence();

  ur2.close();

}

