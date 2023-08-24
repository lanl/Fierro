/** @class */
/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#include "test_fft3d.h"

template<typename backend_tag, typename scalar_type>
using rcontainer = typename heffte::rtransform<backend_tag>::template buffer_container<scalar_type>;

template<typename backend_tag>
void check_cpu_compile_types(){
    // this code will never run, but should compile
    std::unique_ptr<heffte::rtransform<backend_tag>> transform;
    static_assert(std::is_same<decltype(transform->backward(std::declval<rcontainer<backend_tag, float>>())), rcontainer<backend_tag, float>>::value,
                  "r2r backward has incorrect type");
    static_assert(std::is_same<decltype(transform->backward(std::declval<rcontainer<backend_tag, double>>())), rcontainer<backend_tag, double>>::value,
                  "r2r backward has incorrect type");
    static_assert(std::is_same<decltype(transform->backward_real(std::declval<rcontainer<backend_tag, float>>())), rcontainer<backend_tag, float>>::value,
                  "r2r backward has incorrect type");
    static_assert(std::is_same<decltype(transform->backward_real(std::declval<rcontainer<backend_tag, double>>())), rcontainer<backend_tag, double>>::value,
                  "r2r backward has incorrect type");
}

template<typename backend_tag, typename scalar_type>
void test_cosine_transform(MPI_Comm comm){
    using tvector = typename heffte::fft3d<backend_tag>::template buffer_container<scalar_type>; // std::vector or cuda::vector

    int const me = mpi::comm_rank(comm);
    int const num_ranks = mpi::comm_size(comm);
    assert(num_ranks == 1 or num_ranks == 2 or num_ranks == 4);
    current_test<scalar_type, using_mpi, backend_tag> name(std::string("-np ") + std::to_string(num_ranks) + "  test cosine", comm);

    box3d<> const world = {{0, 0, 0}, {1, 2, 3}};
    std::vector<scalar_type> world_input(world.count());
    std::iota(world_input.begin(), world_input.end(), 1.0);
    std::vector<scalar_type> world_result = [&]()->std::vector<scalar_type>{
        if (std::is_same<backend_tag, backend::stock_cos>::value or std::is_same<backend_tag, backend::fftw_cos>::value
            or std::is_same<backend_tag, backend::mkl_cos>::value or std::is_same<backend_tag, backend::cufft_cos>::value
            or std::is_same<backend_tag, backend::rocfft_cos>::value or std::is_same<backend_tag, backend::onemkl_cos>::value) {
            return std::vector<scalar_type>{2.4e+03, -6.7882250993908571e+01, -2.2170250336881628e+02, 0.0, 0.0, 0.0, -9.0844474461089760e+02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -6.4561180200187039e+01, 0.0, 0.0, 0.0, 0.0, 0.0};
        } else if (std::is_same<backend_tag, backend::stock_sin>::value or std::is_same<backend_tag, backend::fftw_sin>::value
            or std::is_same<backend_tag, backend::mkl_sin>::value or std::is_same<backend_tag, backend::cufft_sin>::value
            or std::is_same<backend_tag, backend::rocfft_sin>::value or std::is_same<backend_tag, backend::onemkl_sin>::value) {
            return std::vector<scalar_type>{7.3910362600902943e+02, -4.1810014876044050e+01, -1.0241320258448191e+02, 0.0, 3.6955181300451477e+02, -2.0905007438022025e+01, -3.8400000000000006e+02, 0.0, 0.0, 0.0, -1.9200000000000003e+02, 0.0, 3.0614674589207186e+02, -1.7318275204678301e+01, -4.2420937476555700e+01, 0.0, 1.5307337294603599e+02, -8.6591376023391504e+00, -2.7152900397563417e+02, 0.0, 0.0, 0.0, -1.3576450198781720e+02, 0.0};
        } else if (std::is_same<backend_tag, backend::fftw_cos1>::value) {
            return std::vector<scalar_type>{600.0, -24.0, -48.0, 0.0, 0.0, 0.0, -192.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -48.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        } else if (std::is_same<backend_tag, backend::fftw_sin1>::value) {
            return std::vector<scalar_type>{1.2869458511849268e+03, -5.1477834047397081e+01, -1.7058253619218675e+02, 0.0, 2.2080499998375967e+02, -8.8321999993503919e+00, -6.9064761758659802e+02, 0.0, -5.1623951549498826e-15, 0.0, -1.1849639753652639e+02, 0.0, 3.0380670424097087e+02, -1.2152268169638841e+01, -4.0269074315674175e+01, 0.0, 5.2124989768007239e+01, -2.0849995907202907e+00, -1.6303978624871638e+02, 0.0, -2.1868256803083264e-14, 0.0, -2.7973204907458850e+01, 0.0};
        }
    }();

    std::vector<box3d<>> boxes = [&]()->std::vector<box3d<>>{
            if (num_ranks == 1){
                return heffte::split_world(world, std::array<int, 3>{1, 1, 1});
            }else if (num_ranks == 2){
                return heffte::split_world(world, std::array<int, 3>{2, 1, 1});
            }else{
                return heffte::split_world(world, std::array<int, 3>{1, 2, 2});
            }
        }();
    assert(boxes.size() == static_cast<size_t>(num_ranks));
    auto local_input = input_maker<backend_tag, scalar_type>::select(world, boxes[me], world_input);
    auto reference = get_subbox(world, boxes[me], world_result);
    auto reference_inv = get_subbox(world, boxes[me], world_input);

    for(auto const options : make_all_options<backend_tag>()){
        if (not options.use_pencils) continue;
        heffte::rtransform<backend_tag> trans_cos(boxes[me], boxes[me], comm, options);
        tvector forward(trans_cos.size_outbox());

        trans_cos.forward(local_input.data(), forward.data());
        tassert(approx(forward, reference));

        tvector inverse(trans_cos.size_inbox());
        trans_cos.backward(forward.data(), inverse.data(), heffte::scale::full);
        tassert(approx(inverse, reference_inv, (std::is_same<scalar_type, float>::value) ? 0.001 : 1.0));
    }
}


void perform_tests(MPI_Comm const comm){
    all_tests<> name("cosine transforms");

    check_cpu_compile_types<backend::stock_cos>();
    check_cpu_compile_types<backend::stock_sin>();
    test_cosine_transform<backend::stock_cos, float>(comm);
    test_cosine_transform<backend::stock_cos, double>(comm);
    test_cosine_transform<backend::stock_sin, float>(comm);
    test_cosine_transform<backend::stock_sin, double>(comm);
    #ifdef Heffte_ENABLE_FFTW
    test_cosine_transform<backend::fftw_cos, float>(comm);
    test_cosine_transform<backend::fftw_cos, double>(comm);
    test_cosine_transform<backend::fftw_sin, float>(comm);
    test_cosine_transform<backend::fftw_sin, double>(comm);
    test_cosine_transform<backend::fftw_cos1, float>(comm);
    test_cosine_transform<backend::fftw_cos1, double>(comm);
    test_cosine_transform<backend::fftw_sin1, float>(comm);
    test_cosine_transform<backend::fftw_sin1, double>(comm);
    #endif
    #ifdef Heffte_ENABLE_MKL
    test_cosine_transform<backend::mkl_cos, float>(comm);
    test_cosine_transform<backend::mkl_cos, double>(comm);
    test_cosine_transform<backend::mkl_sin, float>(comm);
    test_cosine_transform<backend::mkl_sin, double>(comm);
    #endif
    #ifdef Heffte_ENABLE_CUDA
    check_cpu_compile_types<backend::cufft_cos>();
    check_cpu_compile_types<backend::cufft_sin>();
    test_cosine_transform<backend::cufft_cos, float>(comm);
    test_cosine_transform<backend::cufft_cos, double>(comm);
    test_cosine_transform<backend::cufft_sin, float>(comm);
    test_cosine_transform<backend::cufft_sin, double>(comm);
    #endif
    #ifdef Heffte_ENABLE_ROCM
    test_cosine_transform<backend::rocfft_cos, float>(comm);
    test_cosine_transform<backend::rocfft_cos, double>(comm);
    test_cosine_transform<backend::rocfft_sin, float>(comm);
    test_cosine_transform<backend::rocfft_sin, double>(comm);
    #endif
    #ifdef Heffte_ENABLE_ONEAPI
    test_cosine_transform<backend::onemkl_cos, float>(comm);
    test_cosine_transform<backend::onemkl_cos, double>(comm);
    test_cosine_transform<backend::onemkl_sin, float>(comm);
    test_cosine_transform<backend::onemkl_sin, double>(comm);
    #endif
}

int main(int argc, char *argv[]){

    MPI_Init(&argc, &argv);

    perform_tests(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
