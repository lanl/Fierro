/** @class */
/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#include "test_fft3d.h"

template<typename backend_tag, bool empty_input_boxes = false>
void test_subcomm_cases(MPI_Comm const comm){
    using location_tag = typename backend::buffer_traits<backend_tag>::location;
    using input_type  = double;
    using cinput_type = std::complex<double>;
    using output_type = std::complex<double>;

    int const me        = mpi::comm_rank(comm);
    int const num_ranks = mpi::comm_size(comm);

    current_test<double, using_mpi, backend_tag> name(std::string("-np ") + std::to_string(num_ranks) + "  test subcommunicator", comm);

    box3d<> const world = {{0, 0, 0}, {23, 23, 23}};

    std::array<int,3> proc_i = heffte::proc_setup_min_surface(world, num_ranks);

    std::vector<box3d<int>> inboxes  = heffte::split_world(world, proc_i);
    box3d<int> inbox = [&]()->box3d<int>{
        if (empty_input_boxes){
            std::vector<box3d<int>> boxes6 = heffte::split_world(world, std::array<int, 3>{2, 3, 1});
            return (me < 6) ? boxes6[me] : box3d<int>(std::array<int, 3>{0, 0, 0}, std::array<int, 3>{-1, -1, -1});
        }else{
            return inboxes[me];
        }
    }();
    box3d<int> outbox = inboxes[me];

    auto world_input = make_data<input_type>(world);
    auto world_fft = forward_fft<backend_tag>(world, world_input);

    auto cworld_input = make_data<output_type>(world);
    auto cworld_fft   = forward_fft<backend_tag>(world, cworld_input);

    auto local_input  = input_maker<backend_tag, input_type>::select(world, inbox, world_input);
    auto clocal_input = input_maker<backend_tag, cinput_type>::select(world, inbox, cworld_input);

    auto local_ref   = get_subbox(world, outbox, world_fft);
    auto clocal_ref  = get_subbox(world, outbox, cworld_fft);

    backend::device_instance<location_tag> device;

    for(auto const &num_subcomm : std::vector<int>{1, 2, 3, 4}){
        for(int variant=0; variant<4; variant++){
            for(auto const &alg : std::vector<reshape_algorithm>{
                reshape_algorithm::alltoall, reshape_algorithm::alltoallv,
                reshape_algorithm::p2p, reshape_algorithm::p2p_plined}){

                heffte::plan_options options = default_options<backend_tag>();

                options.use_num_subranks(num_subcomm);
                options.use_pencils = (variant / 2 == 0);
                options.use_reorder = (variant % 2 == 0);
                options.algorithm = alg;

                auto fft = make_fft3d<backend_tag>(inbox, outbox, comm, options);

                auto lresult = make_buffer_container<output_type>(device.stream(), fft.size_outbox());
                auto lback   = make_buffer_container<input_type>(device.stream(), fft.size_inbox());
                auto clback  = make_buffer_container<output_type>(device.stream(), fft.size_inbox());

                fft.forward(local_input.data(), lresult.data());
                tassert(approx(lresult, local_ref));

                fft.backward(lresult.data(), lback.data(), heffte::scale::full);
                tassert(approx(local_input, lback));

                fft.forward(clocal_input.data(), lresult.data());
                tassert(approx(lresult, clocal_ref));

                fft.backward(lresult.data(), clback.data(), heffte::scale::full);
                tassert(approx(clocal_input, clback));
            }
        }
    }
}
template<typename backend_tag, bool use_communicator>
void test_subcomm_cases_r2c(MPI_Comm const comm){
    using location_tag = typename backend::buffer_traits<backend_tag>::location;
    using input_type  = float;
    using output_type = std::complex<float>;

    int const me        = mpi::comm_rank(comm);
    int const num_ranks = mpi::comm_size(comm);

    current_test<double, using_mpi, backend_tag> name(std::string("-np ") + std::to_string(num_ranks) + "  test subcommunicator", comm);

    box3d<> const world = {{0, 0, 0}, {23, 23, 23}};
    MPI_Comm subcomm = MPI_COMM_NULL;
    if (use_communicator){
        assert(mpi::comm_size(comm) == 8);
        subcomm = mpi::new_comm_from_group({1, 3, 5, 7}, comm);
    }

    std::array<int,3> proc_i = heffte::proc_setup_min_surface(world, num_ranks);

    std::vector<box3d<int>> inboxes  = heffte::split_world(world, proc_i);

    auto world_input = make_data<input_type>(world);
    auto world_fft = forward_fft<backend_tag>(world, world_input);


    auto local_input  = input_maker<backend_tag, input_type>::select(world, inboxes[me], world_input);

    backend::device_instance<location_tag> device;

    for(auto const &num_subcomm : std::vector<int>{1, 2, 4}){
        for(int dim=0; dim<3; dim++){
            std::vector<box3d<>> cboxes = heffte::split_world(world.r2c(dim), proc_i);
            auto world_fft_r2c = get_subbox(world, world.r2c(dim), world_fft);
            auto local_ref = get_subbox(world.r2c(dim), cboxes[me], world_fft_r2c);

        for(int variant=0; variant<4; variant++){
            for(auto const &alg : std::vector<reshape_algorithm>{
                reshape_algorithm::alltoall, reshape_algorithm::alltoallv,
                reshape_algorithm::p2p, reshape_algorithm::p2p_plined}){

                heffte::plan_options options = default_options<backend_tag>();

                if (use_communicator){
                    options.use_num_subranks(num_subcomm);
                }else{
                    options.use_subcomm(subcomm);
                }
                options.use_pencils = (variant / 2 == 0);
                options.use_reorder = (variant % 2 == 0);
                options.algorithm = alg;

                auto fft = make_fft3d_r2c<backend_tag>(inboxes[me], cboxes[me], dim, comm, options);

                auto lresult = make_buffer_container<output_type>(device.stream(), fft.size_outbox());
                auto lback   = make_buffer_container<input_type>(device.stream(), fft.size_inbox());

                fft.forward(local_input.data(), lresult.data());
                tassert(approx(lresult, local_ref));

                fft.backward(lresult.data(), lback.data(), heffte::scale::full);
                tassert(approx(local_input, lback));
            }
        }
        }
    }
    if (subcomm != MPI_COMM_NULL) mpi::comm_free(subcomm);
}

template<typename backend_tag>
void test_all_subcases(MPI_Comm const comm){
    constexpr bool use_empty_input_boxes = true;
    constexpr bool no_empty_input_boxes = false;

    constexpr bool use_subcomm = true;
    constexpr bool use_num_subcomm = false;

    test_subcomm_cases<backend_tag, no_empty_input_boxes>(comm);
    test_subcomm_cases_r2c<backend_tag, use_num_subcomm>(comm);
    if (mpi::comm_size(comm) == 8){
        test_subcomm_cases<backend_tag, use_empty_input_boxes>(comm);
        test_subcomm_cases_r2c<backend_tag, use_subcomm>(comm);
    }
}

void perform_tests(MPI_Comm const comm){
    all_tests<> name("heffte::fft subcommunicators");

    test_all_subcases<backend::stock>(comm);
    #ifdef Heffte_ENABLE_FFTW
    test_all_subcases<backend::fftw>(comm);
    #endif
    #ifdef Heffte_ENABLE_MKL
    test_all_subcases<backend::mkl>(comm);
    #endif
    #ifdef Heffte_ENABLE_CUDA
    test_all_subcases<backend::cufft>(comm);
    #endif
    #ifdef Heffte_ENABLE_ROCM
    test_all_subcases<backend::rocfft>(comm);
    #endif
    #ifdef Heffte_ENABLE_ONEAPI
    test_all_subcases<backend::onemkl>(comm);
    #endif
}

int main(int argc, char *argv[]){

    MPI_Init(&argc, &argv);

    perform_tests(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
