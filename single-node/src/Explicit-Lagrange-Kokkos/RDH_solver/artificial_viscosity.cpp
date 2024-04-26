#include "ref_elem.h"
#include "mesh.h"
#include "state.h"

void get_artificial_viscosity(CArrayKokkos <double> &sigma_a,
                              const DViewCArrayKokkos <double> &vel,
                              const DViewCArrayKokkos <double> &den,
                              const DViewCArrayKokkos <double> &sspd,
                              const DViewCArrayKokkos <double> &vol,
                              const mesh_t &mesh,
                              const fe_ref_elem_t &ref_elem,
                              const size_t stage){


        // vareps = 0.5*( \nabla u + (\nabla u)^T ) //
        CArrayKokkos <double> eps(mesh.num_elems*mesh.num_leg_gauss_in_elem, mesh.num_dims, mesh.num_dims);

        // Tr(\vareps) //
        CArrayKokkos <double> Tr_eps(mesh.num_elems*mesh.num_leg_gauss_in_elem);

        // mu =  \rho \ell ( c_s + \ell |\Delta u| ) //
        CArrayKokkos <double> mu(mesh.num_elems*mesh.num_leg_gauss_in_elem);

        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);

                mu(gauss_gid) = 0.0;
                Tr_eps(gauss_gid) = 0.0;

                for(int i = 0; i < mesh.num_dims; i++){
                    for (int j = 0; j < mesh.num_dims; j++){
                        
                        eps(gauss_gid, i, j) = 0.0;
                    }// j
                }// i

            }// gauss_lid
        });// elem_gid

        // fill eps
        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                for(int i = 0; i < mesh.num_dims; i++){
                    for (int j = 0; j < mesh.num_dims; j++){
                        for (int dof = 0; dof < mesh.num_nodes_in_elem; dof++){
                            int dof_gid = mesh.nodes_in_elem(elem_gid, dof);
                            eps(gauss_gid, i, j) = 0.5*(ref_elem.gauss_leg_grad_basis(gauss_lid, dof, i)*vel(stage, dof_gid, j )
                                                     + ref_elem.gauss_leg_grad_basis(gauss_lid, dof, j)*vel(stage, dof_gid, i ) );
                        }// dof
                    }// j
                }// i
            }// gauss_lid
        });// elem_gid

        // fill Tr_eps
        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                for(int i = 0; i < mesh.num_dims; i++){
                    Tr_eps(gauss_gid) += eps(gauss_gid, i, i);
                }// i
            }// gauss_lid
        });// elem_gid

        // fill mu
        FOR_ALL(elem_gid, 0, mesh.num_elems,{

            double h = pow( vol(elem_gid), 1.0/mesh.num_dims);
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                
                double coeff = 0.0;
                if (h*Tr_eps(gauss_gid) < 0.0){
                    coeff = 0.33;
                }
                mu(gauss_gid) = den(gauss_gid)*h*( coeff*sspd(gauss_gid) + abs(h*Tr_eps(gauss_gid)) );


            }// gauss_lid
        });// elem_gid

        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);
                for(int i = 0; i < mesh.num_dims; i++){
                    for(int j = 0; j < mesh.num_dims; j++){
                        sigma_a(stage, gauss_gid, i, j) = mu(gauss_gid)*eps(gauss_gid, i, j);
                    }
                }
            }
        });

}// end get_artificial_viscosity

void append_artificial_viscosity(DViewCArrayKokkos <double> &sigma,
                                 const CArrayKokkos <double> &sigma_a,
                                 const mesh_t &mesh,
                                 const size_t stage){

        FOR_ALL(elem_gid, 0, mesh.num_elems,{
            for (int gauss_lid = 0; gauss_lid < mesh.num_leg_gauss_in_elem; gauss_lid++){
                int gauss_gid = mesh.legendre_in_elem(elem_gid, gauss_lid);

                for(int i = 0; i < mesh.num_dims; i++){
                    for(int j = 0; j < mesh.num_dims; j++){

                        sigma(stage, gauss_gid, i, j) += sigma_a(stage, gauss_gid, i, j);

                    }// j
                }// i

            }// gauss_lid
        });// elem_gid

}// end append_artificial_viscosity