#ifndef STATE_H
#define STATE_H  

#include "utilities.h"
#include "matar.h"
#include "geometry.h"

using namespace utils;


/// WARNING: NULLIFY ALL POINTERS AND CREATE DECONTSTRUCTOR
// int* a = new int[10]();


class node_t {

private:
    int num_dim_;
    int num_rk_;

    // Number of nodes
    int num_nodes_;

    // c_array_t <real_t> node_field_;

    real_t * node_field_ = NULL;

public:

    void init_node_state (int num_dim, swage::mesh_t& mesh, int num_rk)
    {

        num_dim_   = num_dim;
        num_rk_    = num_rk;
        num_nodes_ = mesh.num_nodes();

        // node_field_ = c_array_t <real_t> (num_nodes_);

        node_field_ = new real_t[num_nodes_]();

    }


    // **** Node State **** //
    // inline real_t& temp(int node_gid) const
    // {
    //     return node_field_(node_gid);
    // }

    inline real_t& field(int node_gid) const
    {
        return node_field_[node_gid];
    }


    // deconstructor WARNING:POSSIBLY REMOVE????
    ~node_t ( ) {
        delete[] node_field_;
    }

};




class mat_pt_t {

private:

    int num_dim_;
    int num_rk_;

    // **** Material Point State **** //
    int num_matpt_;
    
    // the material region id 
    int *mat_id_ = NULL;     

    // physical position of material point
    real_t *coords_ = NULL;

    // Some field
    real_t *field_ = NULL;

    // Volume
    real_t *volume_ = NULL;


public:

    void init_mat_pt_state (int num_dim, swage::mesh_t& mesh, int num_rk)
    {

        num_dim_   = num_dim;
        num_rk_    = num_rk;
        
        // **** Material Point State **** //
        num_matpt_ = mesh.num_cells();
        mat_id_ = new int[num_matpt_]();

        coords_ = new real_t[num_rk_*num_matpt_*num_dim_]();
        
        field_ = new real_t[num_matpt_]();

        volume_ = new real_t[num_matpt_]();

    }

    inline int& mat_id(int mat_pt_gid) const
    {
        return mat_id_[mat_pt_gid];
    }

    inline real_t& coords(int rk_stage, int mat_pt_gid, int this_dim) const
    {
        return coords_[rk_stage*num_matpt_*num_dim_ + mat_pt_gid*num_dim_ + this_dim];
    }

    inline real_t& field(int mat_pt_gid) const
    {
        return field_[mat_pt_gid];
    }

    inline real_t& volume(int mat_pt_gid) const
    {
        return volume_[mat_pt_gid];
    }

    // deconstructor
    ~mat_pt_t ( ) {

        delete[] mat_id_;
        delete[] coords_;
        delete[] field_;
        delete[] volume_;

    }
};






#endif // end STATE_H
