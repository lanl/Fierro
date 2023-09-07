/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#include "heffte_config.h"

#ifdef Heffte_ENABLE_MAGMA

#include "heffte_magma_helpers.h"

#ifdef Heffte_ENABLE_CUDA
#include <cuda_runtime_api.h>
#include <cuda.h>
#else
#ifndef __HIP_PLATFORM_HCC__
#define __HIP_PLATFORM_HCC__
#endif
#include <hip/hip_runtime.h>
#define HAVE_HIP
#endif
#include "magma_v2.h"

namespace heffte {

namespace gpu {

magma_handle<tag::gpu>::magma_handle(void *gpu_stream){
    magma_init();
    int device;
    #ifdef Heffte_ENABLE_CUDA
    cudaStream_t stream = reinterpret_cast<cudaStream_t>(gpu_stream);
    cudaError_t status = cudaGetDevice(&device);
    if (status != cudaSuccess)
        throw std::runtime_error(std::string("cudaGetDevice() failed with message: ") + cudaGetErrorString(status));
    magma_queue_create_from_cuda(device, stream, nullptr, nullptr, reinterpret_cast<magma_queue_t*>(&handle));
    #else
    hipStream_t stream = reinterpret_cast<hipStream_t>(gpu_stream);
    hipError_t status = hipGetDevice(&device);
    if (status != hipSuccess)
        throw std::runtime_error(std::string("hipGetDevice() failed with message: ") + hipGetErrorString(status));
    magma_queue_create_from_hip(device, stream, nullptr, nullptr, reinterpret_cast<magma_queue_t*>(&handle));
    #endif
}
magma_handle<tag::gpu>::~magma_handle(){
    magma_queue_destroy(reinterpret_cast<magma_queue_t>(handle));
    magma_finalize();
}
void magma_handle<tag::gpu>::scal(int num_entries, double scale_factor, float *data) const{
    magma_sscal(num_entries, static_cast<float>(scale_factor), data, 1, reinterpret_cast<magma_queue_t>(handle));
}
void magma_handle<tag::gpu>::scal(int num_entries, double scale_factor, double *data) const{
    magma_dscal(num_entries, scale_factor, data, 1, reinterpret_cast<magma_queue_t>(handle));
}

}

} // namespace heffte

#endif
