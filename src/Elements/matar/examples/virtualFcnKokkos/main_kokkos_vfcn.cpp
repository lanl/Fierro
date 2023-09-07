#include <iostream>
#include "inherited_inits.hpp"

int main() {


    // Kokkos GPU test
   
    Kokkos::initialize();
    {
    

    int num_parent = 2; // number of materials
    Parent1D parent("parent", num_parent); // Initialize Kokkos View on the GPU of type material, size num_parent
    auto h_parent = Kokkos::create_mirror_view(parent); // Create a host view of the Kokkos View


    AllocateHost(h_parent, 0, BABY2_SIZE); // Function performed on Host to do raw Kokkos allocation of baby2 GPU space inside of Host data structure
    AllocateHost(h_parent, 1, BABY1_SIZE); // Function performed on Host to do raw Kokkos allocation of baby1 GPU space inside of Host data structure

    Kokkos::deep_copy(parent, h_parent); // deep copy Host data (allocated above) to the GPU Kokkos View. GPU View now has the class space allocated

    InitChildModels(parent, 0, baby2{}); // Kokkos Function to create new instances of the baby2 model on the GPU
    InitChildModels(parent, 1, baby1{1.4,1.0}); // Kokkos Function to create new instances of the baby1 models on the GPU

    // Model test, also shows a Kokkos reduction
    double value_1;
    Kokkos::parallel_reduce(
        "CheckValues",                             
        num_parent,                   
        KOKKOS_LAMBDA(const int idx, real_t &lsum)
        {
             lsum += parent(idx).child->math(2.0, 4.0);
        }
        , value_1);                                

    printf("value %f\n", value_1);

    ClearDeviceModels(parent); // Kokkos Function to call deconstructors of objects on the GPU

    FreeHost(h_parent); // Function performed on Host to free the allocated GPU classes inside of the Host mirror

   
    }
    Kokkos::finalize();


    printf("--- finished ---\n");

    return 0;
}
