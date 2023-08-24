/** @class */
/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/
#ifndef HEFFTE_TEST_FFT
#define HEFFTE_TEST_FFT

#include <random>

#include "test_common.h"

#ifdef Heffte_ENABLE_CUDA
#include <cuda_runtime_api.h>
#include <cuda.h>
#endif

template<typename scalar_type, typename index>
std::vector<scalar_type> make_data(box3d<index> const world){
    std::minstd_rand park_miller(4242);
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    std::vector<scalar_type> result(world.count());
    for(auto &r : result)
        r = static_cast<scalar_type>(unif(park_miller));
    return result;
}
template<typename scalar_type, typename index>
std::vector<scalar_type> make_data(int num_batch, box3d<index> const world){
    std::minstd_rand park_miller(4242);
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    std::vector<scalar_type> result(num_batch * world.count());
    for(auto &r : result)
        r = static_cast<scalar_type>(unif(park_miller));
    return result;
}

template<typename scalar_type>
void get_subbox_array(box3d<> const world, box3d<> const box, scalar_type const *input, scalar_type *result){
    int const mid_world = world.size[0];
    int const slow_world = world.size[0] * world.size[1];
    int const mid_box = box.size[0];
    int const slow_box = box.size[0] * box.size[1];
    for(int k=box.low[2]; k<=box.high[2]; k++){
        for(int j=box.low[1]; j<=box.high[1]; j++){
            for(int i=box.low[0]; i<=box.high[0]; i++){
                result[(k - box.low[2]) * slow_box + (j - box.low[1]) * mid_box + (i - box.low[0])]
                    = input[k * slow_world + j * mid_world + i];
            }
        }
    }
}

template<typename scalar_type>
std::vector<scalar_type> get_subbox(box3d<> const world, box3d<> const box, std::vector<scalar_type> const &input){
    std::vector<scalar_type> result(box.count());
    get_subbox_array(world, box, input.data(), result.data());
    return result;
}

template<typename scalar_type>
std::vector<scalar_type> get_subboxes(int num_batch, box3d<> const world, box3d<> const box, std::vector<scalar_type> const &input){
    std::vector<scalar_type> result(num_batch * box.count());
    for(int j=0; j<num_batch; j++){
        get_subbox_array(world, box, input.data() + j * world.count(), result.data() + j * box.count());
    }
    return result;
}

template<typename scalar_type>
std::vector<scalar_type> rescale(box3d<> const world, std::vector<scalar_type> const &data, scale scaling){
    std::vector<scalar_type> result = data;
    double scaling_factor = (scaling == scale::none) ? 1.0 : 1.0 / static_cast<double>(world.count());
    if (scaling == scale::symmetric) scaling_factor = std::sqrt(scaling_factor);
    if (scaling != scale::none)
        data_scaling::apply(result.size(), result.data(), scaling_factor);
    return result;
}

template<typename scalar_type>
std::vector<typename fft_output<scalar_type>::type> convert_to_output(std::vector<scalar_type> const &data){
    std::vector<typename fft_output<scalar_type>::type> result(data.size());
    for(size_t i=0; i<data.size(); i++) result[i] = static_cast<typename fft_output<scalar_type>::type>(data[i]);
    return result;
}

template<typename scalar_type>
std::vector<typename fft_output<scalar_type>::type> get_complex_subbox(box3d<> const world, box3d<> const box, std::vector<scalar_type> const &input){
    return get_subbox(world, box, convert_to_output(input));
}

template<typename backend_tag, typename scalar_type, typename = void>
struct input_maker{
    static std::vector<scalar_type> select(box3d<> const world, box3d<> const box, std::vector<scalar_type> const &input){
        return get_subbox(world, box, input);
    }
    static std::vector<scalar_type> select(int num_batch, box3d<> const world, box3d<> const box,
                                           std::vector<scalar_type> const &input){
        return get_subboxes(num_batch, world, box, input);
    }
};

#ifdef Heffte_ENABLE_GPU
template<typename backend_tag, typename scalar_type>
struct input_maker<backend_tag, scalar_type, typename std::enable_if<backend::uses_gpu<backend_tag>::value, void>::type>{
    template<typename index>
    static gpu::vector<scalar_type> select(box3d<index> const world, box3d<index> const box, std::vector<scalar_type> const &input){
        return gpu::transfer().load(get_subbox(world, box, input));
    }
    template<typename index>
    static gpu::vector<scalar_type> select(int num_batch, box3d<index> const world, box3d<index> const box,
                                           std::vector<scalar_type> const &input){
        return gpu::transfer().load(get_subboxes(num_batch, world, box, input));
    }
};
template<typename scalar_type, typename index>
std::vector<scalar_type> rescale(box3d<index> const world, gpu::vector<scalar_type> const &data, scale scaling){
    return rescale(world, gpu::transfer::unload(data), scaling);
}
#endif

template<typename backend_tag, typename precision_type, typename index>
std::vector<std::complex<precision_type>> forward_fft(box3d<index> const world,
                                                      std::vector<precision_type> const &input, int batch_size = 1){
    auto loaded_input = test_traits<backend_tag>::load(input);
    backend::device_instance<typename backend::buffer_traits<backend_tag>::location> device;
    typename test_traits<backend_tag>::template container<std::complex<precision_type>> loaded_result(input.size());
    auto fft0 = make_executor<backend_tag>(device.stream(), world, 0);
    auto fft1 = make_executor<backend_tag>(device.stream(), world, 1);
    auto fft2 = make_executor<backend_tag>(device.stream(), world, 2);
    auto workspace = make_buffer_container<std::complex<precision_type>>(device.stream(),
                        get_max_work_size(std::array<typename one_dim_backend<backend_tag>::executor*, 3>{fft0.get(), fft1.get(), fft2.get()}));
    for(int j=0; j<batch_size; j++){
        fft0->forward(loaded_input.data() + j * world.count(), loaded_result.data() + j * world.count(), workspace.data());
        fft1->forward(loaded_result.data() + j * world.count(), workspace.data());
        fft2->forward(loaded_result.data() + j * world.count(), workspace.data());
    }
    return test_traits<backend_tag>::unload(loaded_result);
}
template<typename backend_tag, typename precision_type, typename index>
std::vector<std::complex<precision_type>> forward_fft(box3d<index> const world,
                                                      std::vector<std::complex<precision_type>> const &input, int batch_size = 1){
    auto loaded_input = test_traits<backend_tag>::load(input);
    backend::device_instance<typename backend::buffer_traits<backend_tag>::location> device;
    for(int i=0; i<3; i++){
        auto fft = make_executor<backend_tag>(device.stream(), world, i);
        auto workspace = make_buffer_container<std::complex<precision_type>>(device.stream(), fft->workspace_size());
        for(int j=0; j<batch_size; j++){
            fft->forward(loaded_input.data() + j * world.count(), workspace.data());
        }
    }
    return test_traits<backend_tag>::unload(loaded_input);
}

template<typename backend_tag>
void test_fft3d_const_dest2(MPI_Comm comm){
    assert(mpi::comm_size(comm) == 2);
    current_test<int, using_mpi, backend_tag> name("constructor heffte::fft3d", comm);
    box3d<> const world = {{0, 0, 0}, {4, 4, 4}};
    std::vector<box3d<>> boxes = heffte::split_world(world, {2, 1, 1});
    int const me = mpi::comm_rank(comm);
    // construct an instance of heffte::fft3d and delete it immediately
    heffte::fft3d<backend_tag> fft(boxes[me], boxes[me], comm);
}

template<typename backend_tag, typename scalar_type, int h0, int h1, int h2>
void test_fft3d_vectors(MPI_Comm comm){
    // works with ranks 6 and 8 only
    int const num_ranks = mpi::comm_size(comm);
    assert(num_ranks == 6 or num_ranks == 8);
    box3d<> const world = {{0, 0, 0}, {h0, h1, h2}};
    std::string fft_dim = (world.is2d()) ? "  test heffte::fft2d" : "  test heffte::fft3d";
    current_test<scalar_type, using_mpi, backend_tag> name(std::string("-np ") + std::to_string(num_ranks) + fft_dim, comm);
    int const me = mpi::comm_rank(comm);
    auto world_input = make_data<scalar_type>(world);
    auto world_complex = convert_to_output(world_input);
    std::vector<decltype(std::real(world_complex[0]))> world_real(world_input.size());
    for(size_t i=0; i<world_real.size(); i++) world_real[i] = std::real(world_input[i]);
    auto world_fft = forward_fft<backend_tag>(world, world_input);

    std::array<heffte::scale, 3> fscale = {heffte::scale::none, heffte::scale::symmetric, heffte::scale::full};
    std::array<heffte::scale, 3> bscale = {heffte::scale::full, heffte::scale::symmetric, heffte::scale::none};

    for(auto const &options : make_all_options2<backend_tag>()){
    for(int i=0; i<3; i++){
        std::array<int, 3> split = {1, 1, 1};
        if (num_ranks == 6){
            split[i] = 2;
            split[(i+1) % 3] = 3;
        }else if (num_ranks == 8){
            split = {2, 2, 2};
        }
        if (world.is2d()){
            assert(num_ranks == 8 and h1 == 0); // 2D test with 8 ranks uses 1 entry for dimension 1, i.e., high-index 1 is 0
            split = {4, 1, 2};
        }
        std::vector<box3d<>> boxes = heffte::split_world(world, split);
        assert(boxes.size() == static_cast<size_t>(num_ranks));

        // get a semi-random inbox and outbox
        // makes sure that the boxes do not have to match
        int iindex = me, oindex = me; // indexes of the input and outboxes
        if (num_ranks == 6){ // shuffle the boxes
            iindex = (me+2) % num_ranks;
            oindex = (me+3) % num_ranks;
        }else if (num_ranks == 8){
            iindex = (me+3) % num_ranks;
            oindex = (me+5) % num_ranks;
        }

        box3d<> const inbox  = boxes[iindex];
        box3d<> const outbox = boxes[oindex];

        auto local_input         = input_maker<backend_tag, scalar_type>::select(world, inbox, world_input);
        auto local_complex_input = get_subbox(world, inbox, world_complex);
        auto local_real_input    = get_subbox(world, inbox, world_real);
        auto reference_fft       = rescale(world, get_subbox(world, outbox, world_fft), fscale[i]);

        heffte::fft3d<backend_tag> fft(inbox, outbox, comm, options);

        auto result = fft.forward(local_input, fscale[i]);
        tassert(approx(result, reference_fft));

        auto backward_cresult = fft.backward(result, bscale[i]);
        auto backward_scaled_cresult = rescale(world, backward_cresult, scale::none);
        tassert(approx(local_complex_input, backward_scaled_cresult));

        auto backward_rresult = fft.backward_real(result, bscale[i]);
        auto backward_scaled_rresult = rescale(world, backward_rresult, scale::none);
        tassert(approx(backward_scaled_rresult, local_real_input));
    }
    } // different option variants
}

template<typename backend_tag, typename scalar_type, int h0, int h1, int h2>
void test_fft3d_queues(MPI_Comm comm){
    // works with ranks 6, same logic as vectors but uses an external device stream/queue
    // used complex numbers only
    int const num_ranks = mpi::comm_size(comm);
    assert(num_ranks == 6);
    box3d<> const world = {{0, 0, 0}, {h0, h1, h2}};
    current_test<scalar_type, using_mpi, backend_tag> name(std::string("-np ") + std::to_string(num_ranks) + "  test heffte::fft3d (stream)", comm);
    int const me = mpi::comm_rank(comm);
    auto world_input = make_data<scalar_type>(world);
    auto world_fft = forward_fft<backend_tag>(world, world_input);

    std::array<heffte::scale, 3> fscale = {heffte::scale::none, heffte::scale::symmetric, heffte::scale::full};
    std::array<heffte::scale, 3> bscale = {heffte::scale::full, heffte::scale::symmetric, heffte::scale::none};

    auto stream = make_stream(backend_tag());

    for(auto const &options : make_all_options<backend_tag>()){
    //if (mpi::world_rank(0)) std::cout << options << std::endl;

    for(int i=0; i<3; i++){
        std::array<int, 3> split = {1, 1, 1};
        if (num_ranks == 6){
            split[i] = 2;
            split[(i+2) % 3] = 3;
        }
        std::vector<box3d<>> boxes = heffte::split_world(world, split);
        // get a semi-random inbox and outbox
        // makes sure that the boxes do not have to match
        int iindex = me, oindex = me; // indexes of the input and outboxes
        if (num_ranks == 6){ // shuffle the boxes
            iindex = (me+1) % num_ranks;
            oindex = (me+3) % num_ranks;
        }

        box3d<> const inbox  = boxes[iindex];
        box3d<> const outbox = boxes[oindex];

        auto local_input         = input_maker<backend_tag, scalar_type>::select(world, inbox, world_input);
        auto reference_fft       = rescale(world, get_subbox(world, outbox, world_fft), fscale[i]);

        sync_stream(stream);
        heffte::fft3d<backend_tag> fft(stream, inbox, outbox, comm, options);
        sync_stream(stream);

        auto result = fft.forward(local_input, fscale[i]);
        tassert(approx(result, reference_fft));

        auto backward_cresult = fft.backward(result, bscale[i]);
        auto backward_scaled_cresult = rescale(world, backward_cresult, scale::none);
        tassert(approx(local_input, backward_scaled_cresult));
    }
    } // different option variants
    free_stream(stream);
}
template<typename backend_tag, typename scalar_type, int h0, int h1, int h2>
void test_fft3d_r2c_queues(MPI_Comm comm){
    // works with ranks 6 and 8 only
    int const num_ranks = mpi::comm_size(comm);
    assert(num_ranks == 6);
    current_test<scalar_type, using_mpi, backend_tag> name(std::string("-np ") + std::to_string(num_ranks) + "  test heffte::fft3d_r2c (stream)", comm);

    double correction = 1.0; // single precision is less stable, especially for larger problems with 12 mpi ranks
    if (std::is_same<scalar_type, float>::value) correction = 1.0E-2;

    int const me = mpi::comm_rank(comm);
    box3d<> const rworld = {{0, 0, 0}, {h0, h1, h2}};
    auto world_input = make_data<scalar_type>(rworld);

    std::array<heffte::scale, 3> fscale = {heffte::scale::none, heffte::scale::symmetric, heffte::scale::full};
    std::array<heffte::scale, 3> bscale = {heffte::scale::full, heffte::scale::symmetric, heffte::scale::none};

    auto stream = make_stream(backend_tag());

    for(auto const &options : make_all_options<backend_tag>()){
    for(int dim = 0; dim < 3; dim++){
        box3d<> const cworld = rworld.r2c(dim);
        auto world_fft     = get_subbox(rworld, cworld, forward_fft<backend_tag>(rworld, world_input));

        for(int i=0; i<3; i++){
            std::array<int, 3> split = {1, 1, 1};
            split[i] = 2;
            split[(i+1) % 3] = 3;

            std::vector<box3d<>> rboxes = heffte::split_world(rworld, split);
            std::vector<box3d<>> cboxes = heffte::split_world(cworld, split);

            assert(rboxes.size() == static_cast<size_t>(num_ranks));
            assert(cboxes.size() == static_cast<size_t>(num_ranks));

            // same as the regular vector API test
            int iindex = me, oindex = me; // indexes of the input and outboxes
            iindex = (me+2) % num_ranks;
            oindex = (me+3) % num_ranks;

            box3d<> const inbox  = rboxes[iindex];
            box3d<> const outbox = cboxes[oindex];

            auto local_input   = input_maker<backend_tag, scalar_type>::select(rworld, inbox, world_input);
            auto reference_fft = rescale(rworld, get_subbox(cworld, outbox, world_fft), fscale[i]);

            sync_stream(stream);
            heffte::fft3d_r2c<backend_tag> fft(stream, inbox, outbox, dim, comm, options);
            sync_stream(stream);

            auto result = fft.forward(local_input, fscale[i]);
            tassert(approx(result, reference_fft, correction));

            auto backward_result = fft.backward(result, bscale[i]);
            auto backward_scaled_result = rescale(rworld, backward_result, scale::none);
            tassert(approx(local_input, backward_scaled_result));
        }
    }
    } // different option variants

    free_stream(stream);
}

template<typename backend_tag, typename scalar_type, int h0, int h1>
void test_fft3d_vectors_2d(MPI_Comm comm){
    // works with ranks 4 and 6 and std::complex<float> or std::complex<double> as scalar_type
    int const num_ranks = mpi::comm_size(comm);
    assert(num_ranks == 4 or num_ranks == 6);
    current_test<scalar_type, using_mpi, backend_tag> name(std::string("-np ") + std::to_string(num_ranks) + "  test heffte::fft2d", comm);
    int const me = mpi::comm_rank(comm);
    box3d<> const world = {{0, 0, 0}, {h0, h1, 0}};
    auto world_input = make_data<scalar_type>(world);
    auto world_fft = forward_fft<backend_tag>(world, world_input);

    std::array<heffte::scale, 3> fscale = {heffte::scale::none, heffte::scale::symmetric, heffte::scale::full};
    std::array<heffte::scale, 3> bscale = {heffte::scale::full, heffte::scale::symmetric, heffte::scale::none};

    for(auto const &options : make_all_options<backend_tag>()){
    for(int i=0; i<3; i++){
        std::array<int, 3> split = (i == 0) ? std::array<int, 3>{3, 2, 1} : std::array<int, 3>{2, 3, 1};
        if (num_ranks == 4) split = {2, 2, 1};
        std::vector<box3d<>> boxes = heffte::split_world(world, split);
        assert(boxes.size() == static_cast<size_t>(num_ranks));

        box2d<> const inbox  = box2d<>(std::array<int, 2>{boxes[me].low[0], boxes[me].low[1]},
                                       std::array<int, 2>{boxes[me].high[0], boxes[me].high[1]});
        box2d<> const outbox = box2d<>(std::array<int, 2>{boxes[me].low[0], boxes[me].low[1]},
                                       std::array<int, 2>{boxes[me].high[0], boxes[me].high[1]});

        auto local_input   = input_maker<backend_tag, scalar_type>::select(world, inbox, world_input);
        auto reference_fft = rescale(world, get_subbox(world, outbox, world_fft), fscale[i]);

        heffte::fft2d<backend_tag> fft(inbox, outbox, comm, options);

        auto result = fft.forward(local_input, fscale[i]);
        tassert(approx(result, reference_fft));

        auto backward_result = fft.backward(result, bscale[i]);
        tassert(approx(backward_result, local_input));
    }
    } // different option variants
}

template<typename backend_tag, typename scalar_type, int h0, int h1, int h2>
void test_fft3d_arrays(MPI_Comm comm){
    using output_type = typename fft_output<scalar_type>::type; // complex type of the output
    using input_container  = typename heffte::fft3d<backend_tag>::template buffer_container<scalar_type>; // std::vector or cuda::vector
    using output_container = typename heffte::fft3d<backend_tag>::template buffer_container<output_type>; // std::vector or cuda::vector

    // works with ranks 2 and 12 only
    int const num_ranks = mpi::comm_size(comm);
    assert(num_ranks == 1 or num_ranks == 2 or num_ranks == 12);
    current_test<scalar_type, using_mpi, backend_tag> name(std::string("-np ") + std::to_string(num_ranks) + "  test heffte::fft3d", comm);

    int const me = mpi::comm_rank(comm);
    box3d<> const world = {{0, 0, 0}, {h0, h1, h2}};
    auto world_input   = make_data<scalar_type>(world);
    auto world_complex = convert_to_output(world_input); // if using real input, convert to the complex output type
    auto world_fft     = forward_fft<backend_tag>(world, world_input); // compute reference fft

    for(auto const &options : make_all_options<backend_tag>()){
    for(int i=0; i<4; i++){ // 3-iterations for different splits, 4-th special one with different in-out boxes for 2 ranks
        // split the world into processors
        std::array<int, 3> split = {1, 1, 1};
        if (num_ranks == 2){
            split[(i < 3) ? i : 2] = 2;
        }else if (num_ranks == 12){
            if (i > 2) continue;
            split = {2, 2, 2};
            split[i] = 3;
        }
        std::vector<box3d<>> boxes = heffte::split_world(world, split);
        assert(boxes.size() == static_cast<size_t>(num_ranks));

        // get the local input as a cuda::vector or std::vector
        auto local_input = input_maker<backend_tag, scalar_type>::select(world, boxes[me], world_input);

        box3d<> outbox = [&]()->box3d<>{
                if (num_ranks == 2 and i == 3){
                    // in the special case, get a different outbox so that algorithm will use only 1 reshape and not shape 0
                    return heffte::split_world(world, std::array<int, 3>{2, 1, 1})[me];
                }else{
                    return boxes[me];
                }
            }();

        heffte::fft3d<backend_tag> fft(boxes[me], outbox, comm, (i < 3) ? options : plan_options(true, reshape_algorithm::alltoallv, false));

        auto reference_fft = get_subbox(world, outbox, world_fft); // reference solution
        output_container forward(fft.size_outbox()); // computed solution
        output_container workspace(fft.size_workspace());

        if (i < 3)
            fft.forward(local_input.data(), forward.data(), workspace.data()); // compute the forward fft
        else
            fft.forward(local_input.data(), forward.data());
        tassert(approx(forward, reference_fft)); // compare to the reference

        #ifdef Heffte_ENABLE_CUDA
        if (std::is_same<backend_tag, backend::cufft>::value){
            scalar_type * uinput;
            output_type *uoutput;
            cudaMallocManaged(reinterpret_cast<void**>(&uinput), local_input.size() * sizeof(scalar_type));
            cudaMallocManaged(reinterpret_cast<void**>(&uoutput), forward.size() * sizeof(output_type));
            cudaMemcpy(uinput, local_input.data(), local_input.size() * sizeof(scalar_type), cudaMemcpyDeviceToDevice);
            fft.forward(uinput, uoutput);
            std::vector<output_type> cpu_result(forward.size());
            cudaMemcpy(cpu_result.data(), uoutput, forward.size() * sizeof(output_type), cudaMemcpyDeviceToHost);
            tassert(approx(forward, cpu_result));
            cudaFree(uinput);
            cudaFree(uoutput);
        }
        #endif

        input_container rbackward(local_input.size()); // compute backward fft using scalar_type
        if (i < 3)
            fft.backward(forward.data(), rbackward.data(), workspace.data());
        else
            fft.backward(forward.data(), rbackward.data());
        auto backward_result = rescale(world, rbackward, scale::full); // always std::vector
        tassert(approx(local_input, backward_result)); // compare with the original input

        output_container cbackward(local_input.size()); // complex backward transform
        fft.backward(forward.data(), cbackward.data(), workspace.data());
        auto cbackward_result = rescale(world, cbackward, scale::full);
        // convert the world to complex numbers and extract the reference sub-box
        tassert(approx(get_complex_subbox(world, boxes[me], world_input),
                       cbackward_result));

        output_container inplace_buffer(std::max(fft.size_inbox(), fft.size_outbox()));
        backend::data_manipulator<typename heffte::fft3d<backend_tag>::location_tag>::copy_n(
            fft.stream(), local_input.data(), fft.size_inbox(), inplace_buffer.data());
        fft.forward(inplace_buffer.data(), inplace_buffer.data());
        output_container inplace_forward(inplace_buffer.data(), inplace_buffer.data() + fft.size_outbox());
        tassert(approx(inplace_forward, reference_fft)); // compare to the reference

        auto inplace_buffer_copy = inplace_buffer;
        fft.backward(inplace_buffer.data(), reinterpret_cast<scalar_type*>(inplace_buffer.data())); // in-place complex-to-real
        rbackward = input_container(reinterpret_cast<scalar_type*>(inplace_buffer.data()),
                                    reinterpret_cast<scalar_type*>(inplace_buffer.data()) + fft.size_inbox());
        backward_result = rescale(world, rbackward, scale::full); // always std::vector
        tassert(approx(local_input, backward_result));

        fft.backward(inplace_buffer_copy.data(), inplace_buffer_copy.data());
        output_container inplace_backward(inplace_buffer_copy.data(), inplace_buffer_copy.data() + fft.size_inbox());
        cbackward_result = rescale(world, inplace_backward, scale::full);
        tassert(approx(get_complex_subbox(world, boxes[me], world_input), cbackward_result));
    }
    } // different option variants
}

template<typename backend_tag>
void test_batch_cases(MPI_Comm const comm){
    using location_tag = typename backend::buffer_traits<backend_tag>::location;
    using input_type  = double;
    using cinput_type = std::complex<double>;
    using output_type = std::complex<double>;

    int const me        = mpi::comm_rank(comm);
    int const num_ranks = mpi::comm_size(comm);
    int const batch_size = 5;

    current_test<double, using_mpi, backend_tag> name(std::string("-np ") + std::to_string(num_ranks) + "  test batch transform", comm);

    box3d<> const  world = {{0, 0, 0}, {7, 8, 9}};

    std::array<int,3> proc_i = heffte::proc_setup_min_surface(world, num_ranks);

    std::vector<box3d<int>> inboxes  = heffte::split_world(world, proc_i);
    box3d<int> inbox  = inboxes[me];
    box3d<int> outbox = inboxes[me];

    auto world_input = make_data<input_type>(batch_size, world);
    auto world_fft = forward_fft<backend_tag>(world, world_input, batch_size);

    auto cworld_input = make_data<output_type>(batch_size, world);
    auto cworld_fft   = forward_fft<backend_tag>(world, cworld_input, batch_size);

    auto local_input  = input_maker<backend_tag, input_type>::select(batch_size, world, inbox, world_input);
    auto clocal_input = input_maker<backend_tag, cinput_type>::select(batch_size, world, inbox, cworld_input);

    auto local_ref   = get_subboxes(batch_size, world, outbox, world_fft);
    auto clocal_ref  = get_subboxes(batch_size, world, outbox, cworld_fft);

    backend::device_instance<location_tag> device;

    for(int variant=0; variant<4; variant++){
        for(auto const &alg : std::vector<reshape_algorithm>{
            reshape_algorithm::alltoall, reshape_algorithm::alltoallv,
            reshape_algorithm::p2p, reshape_algorithm::p2p_plined}){

            heffte::plan_options options = default_options<backend_tag>();

            options.use_pencils = (variant / 2 == 0);
            options.use_reorder = (variant % 2 == 0);
            options.algorithm = alg;

            auto fft = make_fft3d<backend_tag>(inbox, outbox, comm, options);

            auto lresult = make_buffer_container<output_type>(device.stream(), batch_size * fft.size_outbox());
            auto lback   = make_buffer_container<input_type>(device.stream(), batch_size * fft.size_inbox());
            auto clback  = make_buffer_container<output_type>(device.stream(), batch_size * fft.size_inbox());

            auto workspace  = make_buffer_container<output_type>(device.stream(), batch_size * fft.size_workspace());

            fft.forward(batch_size, local_input.data(), lresult.data(), workspace.data());
            tassert(approx(lresult, local_ref));

            fft.backward(batch_size, lresult.data(), lback.data(), workspace.data(), heffte::scale::full);
            tassert(approx(local_input, lback));

            fft.forward(batch_size, clocal_input.data(), lresult.data());
            tassert(approx(lresult, clocal_ref));

            fft.backward(batch_size, lresult.data(), clback.data(), heffte::scale::full);
            tassert(approx(clocal_input, clback));
        }
    }
}

#endif
