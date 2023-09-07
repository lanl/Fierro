/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#include "heffte_backend_cuda.h"

#ifdef Heffte_ENABLE_MAGMA
#include <cublas.h>
#include "magma_v2.h"
#endif

namespace heffte {

namespace gpu {

void device_set(int active_device){
    if (active_device < 0 or active_device > device_count())
        throw std::runtime_error("device_set() called with invalid cuda device id");
    cuda::check_error(cudaSetDevice(active_device), "cudaSetDevice()");
}

void synchronize_default_stream(){
    cuda::check_error(cudaStreamSynchronize(nullptr), "device sync"); // sync the default stream
}

int device_count(){
    int count;
    cuda::check_error(cudaGetDeviceCount(&count), "cudaGetDeviceCount()" );
    return count;
}

}

namespace cuda {

/*
 * Launch with one thread per entry.
 *
 * If to_complex is true, convert one real number from source to two real numbers in destination.
 * If to_complex is false, convert two real numbers from source to one real number in destination.
 */
template<typename scalar_type, int num_threads, bool to_complex, typename index>
__global__ void real_complex_convert(index num_entries, scalar_type const source[], scalar_type destination[]){
    index i = blockIdx.x * num_threads + threadIdx.x;
    while(i < num_entries){
        if (to_complex){
            destination[2*i] = source[i];
            destination[2*i + 1] = 0.0;
        }else{
            destination[i] = source[2*i];
        }
        i += num_threads * gridDim.x;
    }
}

/*
 * Launch this with one block per line.
 */
template<typename scalar_type, int num_threads, int tuple_size, bool pack, typename index>
__global__ void direct_packer(index nfast, index nmid, index nslow, index line_stride, index plane_stide,
                              scalar_type const source[], scalar_type destination[]){
    index block_index = blockIdx.x;
    while(block_index < nmid * nslow){

        index mid = block_index % nmid;
        index slow = block_index / nmid;

        scalar_type const *block_source = (pack) ?
                            &source[tuple_size * (mid * line_stride + slow * plane_stide)] :
                            &source[block_index * nfast * tuple_size];
        scalar_type *block_destination = (pack) ?
                            &destination[block_index * nfast * tuple_size] :
                            &destination[tuple_size * (mid * line_stride + slow * plane_stide)];

        index i = threadIdx.x;
        while(i < nfast * tuple_size){
            block_destination[i] = block_source[i];
            i += num_threads;
        }

        block_index += gridDim.x;
    }
}

/*
 * Launch this with one block per line of the destination.
 */
template<typename scalar_type, int num_threads, int tuple_size, int map0, int map1, int map2, typename index>
__global__ void transpose_unpacker(index nfast, index nmid, index nslow, index line_stride, index plane_stide,
                                   index buff_line_stride, index buff_plane_stride,
                                   scalar_type const source[], scalar_type destination[]){

    index block_index = blockIdx.x;
    while(block_index < nmid * nslow){

        index j = block_index % nmid;
        index k = block_index / nmid;

        index i = threadIdx.x;
        while(i < nfast){
            if (map0 == 0 and map1 == 1 and map2 == 2){
                destination[tuple_size * (k * plane_stide + j * line_stride + i)] = source[tuple_size * (k * buff_plane_stride + j * buff_line_stride + i)];
                if (tuple_size > 1)
                    destination[tuple_size * (k * plane_stide + j * line_stride + i) + 1] = source[tuple_size * (k * buff_plane_stride + j * buff_line_stride + i) + 1];
            }else if (map0 == 0 and map1 == 2 and map2 == 1){
                destination[tuple_size * (k * plane_stide + j * line_stride + i)] = source[tuple_size * (j * buff_plane_stride + k * buff_line_stride + i)];
                if (tuple_size > 1)
                    destination[tuple_size * (k * plane_stide + j * line_stride + i) + 1] = source[tuple_size * (j * buff_plane_stride + k * buff_line_stride + i) + 1];
            }else if (map0 == 1 and map1 == 0 and map2 == 2){
                destination[tuple_size * (k * plane_stide + j * line_stride + i)] = source[tuple_size * (k * buff_plane_stride + i * buff_line_stride + j)];
                if (tuple_size > 1)
                    destination[tuple_size * (k * plane_stide + j * line_stride + i) + 1] = source[tuple_size * (k * buff_plane_stride + i * buff_line_stride + j) + 1];
            }else if (map0 == 1 and map1 == 2 and map2 == 0){
                destination[tuple_size * (k * plane_stide + j * line_stride + i)] = source[tuple_size * (i * buff_plane_stride + k * buff_line_stride + j)];
                if (tuple_size > 1)
                    destination[tuple_size * (k * plane_stide + j * line_stride + i) + 1] = source[tuple_size * (i * buff_plane_stride + k * buff_line_stride + j) + 1];
            }else if (map0 == 2 and map1 == 1 and map2 == 0){
                destination[tuple_size * (k * plane_stide + j * line_stride + i)] = source[tuple_size * (i * buff_plane_stride + j * buff_line_stride + k)];
                if (tuple_size > 1)
                    destination[tuple_size * (k * plane_stide + j * line_stride + i) + 1] = source[tuple_size * (i * buff_plane_stride + j * buff_line_stride + k) + 1];
            }else if (map0 == 2 and map1 == 0 and map2 == 1){
                destination[tuple_size * (k * plane_stide + j * line_stride + i)] = source[tuple_size * (j * buff_plane_stride + i * buff_line_stride + k)];
                if (tuple_size > 1)
                    destination[tuple_size * (k * plane_stide + j * line_stride + i) + 1] = source[tuple_size * (j * buff_plane_stride + i * buff_line_stride + k) + 1];
            }
            i += num_threads;
        }

        block_index += gridDim.x;
    }
}

/*
 * Call with one thread per entry.
 */
template<typename scalar_type, int num_threads, typename index>
__global__ void simple_scal(index num_entries, scalar_type data[], scalar_type scaling_factor){
    index i = blockIdx.x * num_threads + threadIdx.x;
    while(i < num_entries){
        data[i] *= scaling_factor;
        i += num_threads * gridDim.x;
    }
}

#define BLK_X 256

// DCT-II (REDFT10)
// even symmetry; even-indexed elements are 0 (size 2 N + 1)
// (a b c) -> (0 a 0 b 0 c   0   c 0 b 0 a 0)
template<typename scalar_type>
__global__ void cos_pre_forward_kernel(int N, scalar_type const *input, scalar_type *fft_signal){
    int ind = blockIdx.x*BLK_X + threadIdx.x;

    if (ind < N) {
        fft_signal[2*ind]   = 0.0;
        fft_signal[2*ind+1] = input[ind];
    }
    fft_signal[2*N] = 0.0;
    if (ind < N) {
        fft_signal[4*N-2*ind]  = 0;
        fft_signal[4*N-2*ind-1]= input[ind];
    }
}
// (c1 c2 c3 ...) -> (c1.x c2.x c3.x)
template<typename scalar_type>
__global__ void cos_post_forward_kernel(int N, scalar_type const *fft_signal, scalar_type *result){
    int ind = blockIdx.x*BLK_X + threadIdx.x;

    if (ind < N) {
        result[ind] = fft_signal[2*ind];
    }
}
// DCT-III is the inverse of DCT-II (or IDCT; REDFT01)
// odd symmetry; imaginary parts are set to 0
// (a b c) -> (a,0 b,0 c,0 0,0 -c,0 -b,0 -a,0)
template<typename scalar_type>
__global__ void cos_pre_backward_kernel(int N, scalar_type const *input, scalar_type *fft_signal){
    int ind = blockIdx.x*BLK_X + threadIdx.x;

    if (ind < N) {
        fft_signal[2*ind]   = input[ind];
        fft_signal[2*ind+1] = 0.0;
    }
    if (ind == 0){
        fft_signal[2*N]   = 0.0;
        fft_signal[2*N+1] = 0.0;
    }
    if (ind < N) {
        fft_signal[2*(N+ind+1)]   = -input[N-ind-1];
        fft_signal[2*(N+ind+1)+1] = 0.0;
    }
}
// extract the odd elements
// (a b c d e f) -> (b d f)
template<typename scalar_type>
__global__ void cos_post_backward_kernel(int N, scalar_type const *fft_signal, scalar_type *result){
    int ind = blockIdx.x*BLK_X + threadIdx.x;

    if (ind < N) {
        result[ind] =  fft_signal[1+2*ind];
    }
}

// DST-II (RODFT10)
// odd symmetry; even-indexed elements are 0 (size 2 N + 1)
// (a b c) -> (0 a 0 b 0 c   0   -c 0 -b 0 -a 0)
template<typename scalar_type>
__global__ void sin_pre_forward_kernel(int N, scalar_type const *input, scalar_type *fft_signal){
    int ind = blockIdx.x*BLK_X + threadIdx.x;

    if (ind < N) {
        fft_signal[2*ind]   = 0;
        fft_signal[2*ind+1] = input[ind];
    }
    fft_signal[2*N] = 0.;
    if (ind < N) {
        fft_signal[4*N-2*ind]  = 0;
        fft_signal[4*N-2*ind-1]= -input[ind];
    }
}

// (c1 c2 c3 ...) -> (-c1.y -c2.y -c3.y)
template<typename scalar_type>
__global__ void sin_post_forward_kernel(int N, scalar_type const *fft_signal, scalar_type *result){
    int ind = blockIdx.x*BLK_X + threadIdx.x;

    if (ind < N)
        result[ind] =  -fft_signal[2*(ind+1)+1];
}

// DST-III is the inverse of DST-II (or IDST; RODFT01)
// even symmetry; real parts are set to 0; size 2N+1
// (a b c) -> (0,0 0,-a 0,-b 0,-c 0,-b 0,-a 0,0)
template<typename scalar_type>
__global__ void sin_pre_backward_kernel(int N, scalar_type const *input, scalar_type *fft_signal){
    int ind = blockIdx.x*BLK_X + threadIdx.x;

    if (ind == 0){
        fft_signal[0] = 0.0;
        fft_signal[1] = 0.0;
    }
    if (ind < N) {
        fft_signal[2*(ind+1)]   = 0.0;
        fft_signal[2*(ind+1)+1] = -input[ind];
    }
    if (ind == N-1){
        fft_signal[4*N]   = 0.0;
        fft_signal[4*N+1] = 0.0;
    } else if (ind < N) {
        fft_signal[2*(N+ind+1)]   = 0.0;
        fft_signal[2*(N+ind+1)+1] = -input[N-ind-2];
    }
}

// extract the odd elements
// (a b c d e f) -> (b d f)
template<typename scalar_type>
__global__ void sin_post_backward_kernel(int N, scalar_type const *fft_signal, scalar_type *result){
    int ind = blockIdx.x*BLK_X + threadIdx.x;

    if (ind < N) {
        result[ind] =  fft_signal[1+2*ind];
    }
}


/*
 * Create a 1-D CUDA thread grid using the total_threads and number of threads per block.
 * Basically, computes the number of blocks but no more than 65536.
 */
struct thread_grid_1d{
    // Compute the threads and blocks.
    thread_grid_1d(int total_threads, int num_per_block) :
        threads(num_per_block),
        blocks(std::min(total_threads / threads + ((total_threads % threads == 0) ? 0 : 1), 65536))
    {}
    // number of threads
    int const threads;
    // number of blocks
    int const blocks;
};

// max number of cuda threads (Volta supports more, but I don't think it matters)
constexpr int max_threads  = 1024;
// allows expressive calls to_complex or not to_complex
constexpr bool to_complex  = true;
// allows expressive calls to_pack or not to_pack
constexpr bool to_pack     = true;

template<typename precision_type, typename index>
void convert(cudaStream_t stream, index num_entries, precision_type const source[], std::complex<precision_type> destination[]){
    thread_grid_1d grid(num_entries, max_threads);
    real_complex_convert<precision_type, max_threads, to_complex><<<grid.blocks, grid.threads, 0, stream>>>(num_entries, source, reinterpret_cast<precision_type*>(destination));
}
template<typename precision_type, typename index>
void convert(cudaStream_t stream, index num_entries, std::complex<precision_type> const source[], precision_type destination[]){
    thread_grid_1d grid(num_entries, max_threads);
    real_complex_convert<precision_type, max_threads, not to_complex><<<grid.blocks, grid.threads, 0, stream>>>(num_entries, reinterpret_cast<precision_type const*>(source), destination);
}

#define heffte_instantiate_convert(precision, index) \
    template void convert<precision, index>(cudaStream_t, index num_entries, precision const source[], std::complex<precision> destination[]); \
    template void convert<precision, index>(cudaStream_t, index num_entries, std::complex<precision> const source[], precision destination[]); \

heffte_instantiate_convert(float, int)
heffte_instantiate_convert(double, int)
heffte_instantiate_convert(float, long long)
heffte_instantiate_convert(double, long long)

/*
 * For float and double, defines type = <float/double> and tuple_size = 1
 * For complex float/double, defines type <float/double> and typle_size = 2
 */
template<typename scalar_type> struct precision{
    using type = scalar_type;
    static const int tuple_size = 1;
};
template<typename precision_type> struct precision<std::complex<precision_type>>{
    using type = precision_type;
    static const int tuple_size = 2;
};

template<typename scalar_type, typename index>
void direct_pack(cudaStream_t stream, index nfast, index nmid, index nslow, index line_stride, index plane_stide,
                 scalar_type const source[], scalar_type destination[]){
    constexpr index max_blocks = 65536;
    using prec = typename precision<scalar_type>::type;
    direct_packer<prec, max_threads, precision<scalar_type>::tuple_size, to_pack>
            <<<std::min(nmid * nslow, max_blocks), max_threads, 0, stream>>>(nfast, nmid, nslow, line_stride, plane_stide,
            reinterpret_cast<prec const*>(source), reinterpret_cast<prec*>(destination));
}

template<typename scalar_type, typename index>
void direct_unpack(cudaStream_t stream, index nfast, index nmid, index nslow, index line_stride, index plane_stide, scalar_type const source[], scalar_type destination[]){
    constexpr index max_blocks = 65536;
    using prec = typename precision<scalar_type>::type;
    direct_packer<prec, max_threads, precision<scalar_type>::tuple_size, not to_pack>
            <<<std::min(nmid * nslow, max_blocks), max_threads, 0, stream>>>(nfast, nmid, nslow, line_stride, plane_stide,
            reinterpret_cast<prec const*>(source), reinterpret_cast<prec*>(destination));
}

template<typename scalar_type, typename index>
void transpose_unpack(cudaStream_t stream, index nfast, index nmid, index nslow, index line_stride, index plane_stride,
                      index buff_line_stride, index buff_plane_stride, int map0, int map1, int map2,
                      scalar_type const source[], scalar_type destination[]){
    constexpr index max_blocks = 65536;
    using prec = typename precision<scalar_type>::type;
    if (map0 == 0 and map1 == 1 and map2 == 2){
        transpose_unpacker<prec, max_threads, precision<scalar_type>::tuple_size, 0, 1, 2>
                <<<std::min(nmid * nslow, max_blocks), max_threads, 0, stream>>>
                (nfast, nmid, nslow, line_stride, plane_stride, buff_line_stride, buff_plane_stride,
                 reinterpret_cast<prec const*>(source), reinterpret_cast<prec*>(destination));
    }else if (map0 == 0 and map1 == 2 and map2 == 1){
        transpose_unpacker<prec, max_threads, precision<scalar_type>::tuple_size, 0, 2, 1>
                <<<std::min(nmid * nslow, max_blocks), max_threads, 0, stream>>>
                (nfast, nmid, nslow, line_stride, plane_stride, buff_line_stride, buff_plane_stride,
                 reinterpret_cast<prec const*>(source), reinterpret_cast<prec*>(destination));
    }else if (map0 == 1 and map1 == 0 and map2 == 2){
        transpose_unpacker<prec, max_threads, precision<scalar_type>::tuple_size, 1, 0, 2>
                <<<std::min(nmid * nslow, max_blocks), max_threads, 0, stream>>>
                (nfast, nmid, nslow, line_stride, plane_stride, buff_line_stride, buff_plane_stride,
                 reinterpret_cast<prec const*>(source), reinterpret_cast<prec*>(destination));
    }else if (map0 == 1 and map1 == 2 and map2 == 0){
        transpose_unpacker<prec, max_threads, precision<scalar_type>::tuple_size, 1, 2, 0>
                <<<std::min(nmid * nslow, max_blocks), max_threads, 0, stream>>>
                (nfast, nmid, nslow, line_stride, plane_stride, buff_line_stride, buff_plane_stride,
                 reinterpret_cast<prec const*>(source), reinterpret_cast<prec*>(destination));
    }else if (map0 == 2 and map1 == 0 and map2 == 1){
        transpose_unpacker<prec, max_threads, precision<scalar_type>::tuple_size, 2, 0, 1>
                <<<std::min(nmid * nslow, max_blocks), max_threads, 0, stream>>>
                (nfast, nmid, nslow, line_stride, plane_stride, buff_line_stride, buff_plane_stride,
                 reinterpret_cast<prec const*>(source), reinterpret_cast<prec*>(destination));
    }else if (map0 == 2 and map1 == 1 and map2 == 0){
        transpose_unpacker<prec, max_threads, precision<scalar_type>::tuple_size, 2, 1, 0>
                <<<std::min(nmid * nslow, max_blocks), max_threads, 0, stream>>>
                (nfast, nmid, nslow, line_stride, plane_stride, buff_line_stride, buff_plane_stride,
                 reinterpret_cast<prec const*>(source), reinterpret_cast<prec*>(destination));
    }
}

#define heffte_instantiate_packers(index) \
template void direct_pack<float, index>(cudaStream_t, index, index, index, index, index, float const source[], float destination[]); \
template void direct_pack<double, index>(cudaStream_t, index, index, index, index, index, double const source[], double destination[]); \
template void direct_pack<std::complex<float>, index>(cudaStream_t, index, index, index, index, index, \
                                                      std::complex<float> const source[], std::complex<float> destination[]); \
template void direct_pack<std::complex<double>, index>(cudaStream_t, index, index, index, index, index, \
                                                       std::complex<double> const source[], std::complex<double> destination[]); \
\
template void direct_unpack<float, index>(cudaStream_t, index, index, index, index, index, float const source[], float destination[]); \
template void direct_unpack<double, index>(cudaStream_t, index, index, index, index, index, double const source[], double destination[]); \
template void direct_unpack<std::complex<float>, index>(cudaStream_t, index, index, index, index, index, \
                                                        std::complex<float> const source[], std::complex<float> destination[]); \
template void direct_unpack<std::complex<double>, index>(cudaStream_t, index, index, index, index, index, \
                                                         std::complex<double> const source[], std::complex<double> destination[]); \
\
template void transpose_unpack<float, index>(cudaStream_t, index, index, index, index, index, index, index, int, int, int, \
                                             float const source[], float destination[]); \
template void transpose_unpack<double, index>(cudaStream_t, index, index, index, index, index, index, index, int, int, int, \
                                              double const source[], double destination[]); \
template void transpose_unpack<std::complex<float>, index>(cudaStream_t, index, index, index, index, index, index, index, int, int, int, \
                                                           std::complex<float> const source[], std::complex<float> destination[]); \
template void transpose_unpack<std::complex<double>, index>(cudaStream_t, index, index, index, index, index, index, index, int, int, int, \
                                                            std::complex<double> const source[], std::complex<double> destination[]); \

heffte_instantiate_packers(int)
heffte_instantiate_packers(long long)


template<typename scalar_type, typename index>
void scale_data(cudaStream_t stream, index num_entries, scalar_type *data, double scale_factor){
    thread_grid_1d grid(num_entries, max_threads);
    simple_scal<scalar_type, max_threads><<<grid.blocks, grid.threads, 0, stream>>>(num_entries, data, static_cast<scalar_type>(scale_factor));
}

template void scale_data<float, int>(cudaStream_t, int num_entries, float *data, double scale_factor);
template void scale_data<double, int>(cudaStream_t, int num_entries, double *data, double scale_factor);
template void scale_data<float, long long>(cudaStream_t, long long num_entries, float *data, double scale_factor);
template void scale_data<double, long long>(cudaStream_t, long long num_entries, double *data, double scale_factor);

template<typename precision>
void cos_pre_pos_processor::pre_forward(cudaStream_t stream, int length, precision const input[], precision fft_signal[]){
    dim3 threads( BLK_X, 1 );
    dim3 grid( (length + BLK_X-1)/BLK_X, 1 );
    cos_pre_forward_kernel<<<grid, threads, 0, stream>>>(length, input, fft_signal);
}
template<typename precision>
void cos_pre_pos_processor::post_forward(cudaStream_t stream, int length, std::complex<precision> const fft_result[], precision result[]){
    dim3 threads( BLK_X, 1 );
    dim3 grid( (length + BLK_X-1)/BLK_X, 1 );
    cos_post_forward_kernel<<<grid, threads, 0, stream>>>(length, reinterpret_cast<precision const*>(fft_result), result);
}
template<typename precision>
void cos_pre_pos_processor::pre_backward(cudaStream_t stream, int length, precision const input[], std::complex<precision> fft_signal[]){
    dim3 threads( BLK_X, 1 );
    dim3 grid( (length + BLK_X-1)/BLK_X, 1 );
    cos_pre_backward_kernel<<<grid, threads, 0, stream>>>(length, input, reinterpret_cast<precision*>(fft_signal));
}
template<typename precision>
void cos_pre_pos_processor::post_backward(cudaStream_t stream, int length, precision const fft_result[], precision result[]){
    dim3 threads( BLK_X, 1 );
    dim3 grid( (length + BLK_X-1)/BLK_X, 1 );
    cos_post_backward_kernel<<<grid, threads, 0, stream>>>(length, fft_result, result);
}
template<typename precision>
void sin_pre_pos_processor::pre_forward(cudaStream_t stream, int length, precision const input[], precision fft_signal[]){
    dim3 threads( BLK_X, 1 );
    dim3 grid( (length + BLK_X-1)/BLK_X, 1 );
    sin_pre_forward_kernel<<<grid, threads, 0, stream>>>(length, input, fft_signal);
}
template<typename precision>
void sin_pre_pos_processor::post_forward(cudaStream_t stream, int length, std::complex<precision> const fft_result[], precision result[]){
    dim3 threads( BLK_X, 1 );
    dim3 grid( (length + BLK_X-1)/BLK_X, 1 );
    sin_post_forward_kernel<<<grid, threads, 0, stream>>>(length, reinterpret_cast<precision const*>(fft_result), result);
}
template<typename precision>
void sin_pre_pos_processor::pre_backward(cudaStream_t stream, int length, precision const input[], std::complex<precision> fft_signal[]){
    dim3 threads( BLK_X, 1 );
    dim3 grid( (length + BLK_X-1)/BLK_X, 1 );
    sin_pre_backward_kernel<<<grid, threads, 0, stream>>>(length, input, reinterpret_cast<precision*>(fft_signal));
}
template<typename precision>
void sin_pre_pos_processor::post_backward(cudaStream_t stream, int length, precision const fft_result[], precision result[]){
    dim3 threads( BLK_X, 1 );
    dim3 grid( (length + BLK_X-1)/BLK_X, 1 );
    sin_post_backward_kernel<<<grid, threads, 0, stream>>>(length, fft_result, result);
}

#define heffte_instantiate_cos(precision) \
    template void cos_pre_pos_processor::pre_forward<precision>(cudaStream_t, int, precision const[], precision[]); \
    template void cos_pre_pos_processor::post_forward<precision>(cudaStream_t, int,  std::complex<precision> const[], precision[]); \
    template void cos_pre_pos_processor::pre_backward<precision>(cudaStream_t, int, precision const[], std::complex<precision>[]); \
    template void cos_pre_pos_processor::post_backward<precision>(cudaStream_t, int, precision const[], precision[]); \
    template void sin_pre_pos_processor::pre_forward<precision>(cudaStream_t, int, precision const[], precision[]); \
    template void sin_pre_pos_processor::post_forward<precision>(cudaStream_t, int,  std::complex<precision> const[], precision[]); \
    template void sin_pre_pos_processor::pre_backward<precision>(cudaStream_t, int, precision const[], std::complex<precision>[]); \
    template void sin_pre_pos_processor::post_backward<precision>(cudaStream_t, int, precision const[], precision[]); \

heffte_instantiate_cos(float)
heffte_instantiate_cos(double)

} // namespace cuda


} // namespace heffte
