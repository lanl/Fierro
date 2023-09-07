/** @class */
/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#include "test_fft3d.h"

template<typename backend_tag, typename scalar_type, int r2c, int h0, int h1, int h2>
void test_int64(MPI_Comm comm){
    // works with ranks 4 and 6 and std::complex<float> or std::complex<double> as scalar_type
    int const num_ranks = mpi::comm_size(comm);
    assert(num_ranks == 4);
    current_test<scalar_type, using_mpi, backend_tag> name(std::string("-np ") + std::to_string(num_ranks) + "  test int/long long", comm);
    int const me = mpi::comm_rank(comm);
    box3d<int> const world = {{0, 0, 0}, {h0, h1, h2}};
    box3d<long long> const world64 = {{0, 0, 0}, {h0, h1, h2}};
    auto world_input = make_data<scalar_type>(world);

    std::array<heffte::scale, 3> fscale = {heffte::scale::none, heffte::scale::symmetric, heffte::scale::full};
    std::array<heffte::scale, 3> bscale = {heffte::scale::full, heffte::scale::symmetric, heffte::scale::none};

    for(auto const &options : make_all_options<backend_tag>()){
    for(int i=0; i<3; i++){
        std::array<int, 3> split = {2, 2, 2};
        split[i] = 1;
        std::vector<box3d<int>> boxes = heffte::split_world(world, split);
        std::vector<box3d<long long>> boxes64 = heffte::split_world(world64, split);
        assert(boxes.size() == static_cast<size_t>(num_ranks));
        auto local_input   = input_maker<backend_tag, scalar_type>::select(world, boxes[me], world_input);

        auto fft = heffte::make_fft3d<backend_tag>(boxes[me], boxes[me], comm, options);
        auto fft64 = heffte::make_fft3d<backend_tag>(boxes64[me], boxes64[me], comm, options);
        static_assert(std::is_same<decltype(fft), heffte::fft3d<backend_tag, int>>::value,
                      "heffte::make_fft3d failed to auto-detect int");
        static_assert(std::is_same<decltype(fft64), heffte::fft3d<backend_tag, long long>>::value,
                      "heffte::make_fft3d failed to auto-detect long long");

        auto result   = fft.forward(local_input, fscale[i]);
        auto result64 = fft64.forward(local_input, fscale[i]);
        tassert(approx(result, result64));

        auto back   = fft.backward_real(result, bscale[i]);
        auto back64 = fft64.backward_real(result, bscale[i]);
        tassert(approx(back, back64));

        std::vector<box3d<int>> boxes_r2c = heffte::split_world(world.r2c(r2c), split);
        std::vector<box3d<long long>> boxes64_r2c = heffte::split_world(world64.r2c(r2c), split);
        auto fft_r2c = heffte::make_fft3d_r2c<backend_tag>(boxes[me], boxes_r2c[me], r2c, comm, options);
        auto fft64_r2c = heffte::make_fft3d_r2c<backend_tag>(boxes64[me], boxes64_r2c[me], r2c, comm, options);
        static_assert(std::is_same<decltype(fft_r2c), heffte::fft3d_r2c<backend_tag, int>>::value,
                      "heffte::make_fft3d_r2c failed to auto-detect int");
        static_assert(std::is_same<decltype(fft64_r2c), heffte::fft3d_r2c<backend_tag, long long>>::value,
                      "heffte::make_fft3d_r2c failed to auto-detect long long");

        result   = fft_r2c.forward(local_input, fscale[i]);
        result64 = fft64_r2c.forward(local_input, fscale[i]);
        tassert(approx(result, result64));

        back   = fft_r2c.backward(result, bscale[i]);
        back64 = fft64_r2c.backward(result, bscale[i]);
        tassert(approx(back, back64));
    }
    } // different option variants
}

template<typename backend_tag>
void test_fft3d_cases(MPI_Comm const comm){
    test_int64<backend_tag, float, 0, 31, 31, 30>(comm);
    test_int64<backend_tag, double, 0, 22, 23, 25>(comm);
    test_int64<backend_tag, float, 1, 11, 21, 14>(comm);
    test_int64<backend_tag, double, 1, 15, 15, 15>(comm);
    test_int64<backend_tag, float, 2, 31, 31, 31>(comm);
    test_int64<backend_tag, double, 2, 20, 20, 20>(comm);
}

void perform_tests(MPI_Comm const comm){
    all_tests<> name("heffte::fft class");

    test_fft3d_cases<backend::stock>(comm);
    #ifdef Heffte_ENABLE_FFTW
    test_fft3d_cases<backend::fftw>(comm);
    #endif
    #ifdef Heffte_ENABLE_MKL
    test_fft3d_cases<backend::mkl>(comm);
    #endif
    #ifdef Heffte_ENABLE_CUDA
    test_fft3d_cases<backend::cufft>(comm);
    #endif
    #ifdef Heffte_ENABLE_ROCM
    test_fft3d_cases<backend::rocfft>(comm);
    #endif
    #ifdef Heffte_ENABLE_ONEAPI
    test_fft3d_cases<backend::onemkl>(comm);
    #endif
}

int main(int argc, char *argv[]){

    MPI_Init(&argc, &argv);

    perform_tests(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
