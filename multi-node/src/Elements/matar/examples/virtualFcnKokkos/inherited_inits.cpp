#include "inherited_inits.hpp"

void AllocateHost(ParentHost1D h_parent, u_int idx, size_t size) 
{

    h_parent(idx).child = (child_models *) kmalloc(size);

}

void FreeHost(ParentHost1D h_parent) 
{

    for (int mem = 0; mem < h_parent.extent(0); mem++) {
        kfree(h_parent(mem).child);
    }
    
}

void InitChildModels(Parent1D parent, u_int idx, baby2 baby2_inp) 
{

    Kokkos::parallel_for(
            "CreateObjects", 1, KOKKOS_LAMBDA(const int&) {
                //CreateChildObjects(parent, baby2_inp, idx);
                new ((baby2 *)parent(idx).child) baby2{baby2_inp};
            });

}


void InitChildModels(Parent1D parent, u_int idx, baby1 baby1_inp) 
{

    Kokkos::parallel_for(
            "CreateObjects", 1, KOKKOS_LAMBDA(const int&) {
                //CreateChildObjects(parent, baby1_inp, idx);
                new ((baby1 *)parent(idx).child) baby1{baby1_inp};
            });

}


void ClearDeviceModels(Parent1D parent)
{

    Kokkos::parallel_for(
            "DestroyObjects", 1, KOKKOS_LAMBDA(const int&) {
              parent(0).child->~child_models();
              parent(1).child->~child_models();
            });

}

