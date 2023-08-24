#ifndef KOKKOS_SIMPLE_H
#define KOKKOS_SIMPLE_H

#include <stdio.h>
#include "kokkos_alias.h"

#define BABY1_SIZE ( sizeof(baby1) )
#define BABY2_SIZE ( sizeof(baby2) )


void AllocateHost(ParentHost1D h_parent, u_int idx, size_t size);

void FreeHost(ParentHost1D h_parent); 

void InitChildModels(Parent1D parent, u_int idx, baby1 baby1_inp); 

void InitChildModels(Parent1D parent, u_int idx, baby2 baby2_inp); 

void ClearDeviceModels(Parent1D parent); // would need modification from user by calling the correct deconstructor 

#endif
