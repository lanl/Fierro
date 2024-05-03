#pragma once

#include "matar.h"
#include "mesh.h"
#include "_debug_tools.h"

using namespace mtr;

struct contact_patch_t
{
    size_t gid;  // global patch id
    CArrayKokkos<size_t> nodes_gid;  // global node ids

    void initialize();

};

struct contact_patches_t
{
    CArrayKokkos<contact_patch_t> contact_patches;  // patches that will be checked for contact
    CArrayKokkos<size_t> patches_gid;  // global patch ids
    size_t num_contact_patches;  // total number of patches that will be checked for contact

    void initialize(const mesh_t &mesh, const CArrayKokkos<size_t> &bdy_contact_patches);

};
