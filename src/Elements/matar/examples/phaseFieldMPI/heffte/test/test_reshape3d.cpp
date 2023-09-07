/** @class */
/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#include "test_fft3d.h"

#ifdef Heffte_ENABLE_FFTW
using default_cpu_backend = heffte::backend::fftw;
#elif defined(Heffte_ENABLE_MKL)
using default_cpu_backend = heffte::backend::mkl;
#else
using default_cpu_backend = heffte::backend::stock;
#endif


/*
 * Simple unit test that checks the operation that gathers boxes across an mpi comm.
 */
void test_boxes(MPI_Comm const comm){
    current_test<> test("heffte::mpi::gather_boxes", comm);

    int const me = mpi::comm_rank(comm);

    std::vector<box3d<>> reference_inboxes;
    std::vector<box3d<>> reference_outboxes;

    for(int i=0; i<mpi::comm_size(comm); i++){
        reference_inboxes.push_back({{i, i+1, i+2}, {i+3, i+4, i+5}});
        reference_outboxes.push_back({{i, i+3, i+5}, {i+7, i+6, i+9}});
    }

    ioboxes<> boxes = mpi::gather_boxes(reference_inboxes[me], reference_outboxes[me], comm);

    tassert(match(boxes.in,  reference_inboxes));
    tassert(match(boxes.out, reference_outboxes));
}

/*
 * Returns a vector of data corresponding to a sub-box of the original world.
 * The entries are floating point numbers (real or complex) but have integer values
 * corresponding to the indexes in the world box.
 * Thus, by checking the indexes, it is easy to check if data was moved correctly
 * from one sub-box to another.
 */
template<typename scalar_type>
std::vector<scalar_type> get_subdata(box3d<> const world, box3d<> const subbox){
    // the entries in the master box go 0, 1, 2, 3, 4 ...
    int const wmidstride  = world.size[0];
    int const wslowstride = world.size[0] * world.size[1];
    int const smidstride  = subbox.size[0];
    int const sslowstride = subbox.size[0] * subbox.size[1];

    std::vector<scalar_type> result(subbox.count());

    for(int k = 0; k < subbox.size[2]; k++){
        for(int j = 0; j < subbox.size[1]; j++){
            for(int i = 0; i < subbox.size[0]; i++){
                result[k * sslowstride + j * smidstride + i]
                    = static_cast<scalar_type>((k + world.low[2] + subbox.low[2]) * wslowstride
                                                + (j + world.low[1] + subbox.low[1]) * wmidstride
                                                + i + world.low[0] + subbox.low[0]);
            }
        }
    }
    return result;
}

template<typename location_tag, reshape_algorithm variant, typename index = int>
std::unique_ptr<reshape3d_base<index>>
make_test_reshape3d(typename backend::device_instance<location_tag>::stream_type q, std::vector<box3d<>> const &input_boxes, std::vector<box3d<>> const &output_boxes, MPI_Comm const comm){
    if (variant == reshape_algorithm::alltoallv){
        return make_reshape3d_alltoallv<location_tag>(q, input_boxes, output_boxes, true, comm);
    }else if (variant == reshape_algorithm::alltoall){
        return make_reshape3d_alltoall<location_tag>(q, input_boxes, output_boxes, true, comm);
    }else if (variant == reshape_algorithm::p2p){
        return make_reshape3d_pointtopoint<location_tag>(q, input_boxes, output_boxes, reshape_algorithm::p2p, true, comm);
    }else if (variant == reshape_algorithm::p2p_plined){
        return make_reshape3d_pointtopoint<location_tag>(q, input_boxes, output_boxes, reshape_algorithm::p2p_plined, true, comm);
    }
}

// splits the world box into a set of boxes with gird given by proc_grid
template<int hfast, int hmid, int hslow, int pfast, int pmid, int pslow, typename scalar_type, typename location_tag, reshape_algorithm variant>
void test_cpu(MPI_Comm const comm){
    /*
     * simple test, create a world of indexes going all the way to hfast, hmid and hslow
     * then split the world into boxes numbering pfast, pmid, and pslow, assume that's what each rank owns
     * then create a new world of pencils and assigns a pencil to each rank (see the shuffle comment)
     * more the data from the original configuration to the new and check against reference data
     */
    current_test<scalar_type, using_mpi, location_tag> test("-np " + std::to_string(mpi::comm_size(comm)) + "  "
                                                           + get_description<variant>(), comm);
    tassert( pfast * pmid * pslow == heffte::mpi::comm_size(comm) );

    int const me = heffte::mpi::comm_rank(comm);
    int const shift = 3;

    box3d<> world = {{0, 0, 0}, {hfast, hmid, hslow}};

    auto boxes   = split_world(world, {pfast, pmid, pslow});
    auto pencils = split_world(world, {pfast,    1, pmid * pslow});

    std::vector<box3d<>> rotate_boxes;
    if (std::is_same<scalar_type, std::complex<float>>::value){
        // shuffle the pencil boxes in some tests to check the case when there is no overlap between inbox and outbox
        // for the 2 by 2 grid, this shuffle ensures no overlap
        for(size_t i=0; i<boxes.size(); i++) rotate_boxes.push_back( pencils[(i + shift) % boxes.size()] );
    }else{
        for(auto b : pencils) rotate_boxes.push_back(b);
    }

    // create caches for a reshape algorithm, including creating a new mpi comm
    auto reshape = make_test_reshape3d<location_tag, variant>(nullptr, boxes, rotate_boxes, comm);
    std::vector<scalar_type> workspace(reshape->size_workspace());

    auto input_data     = get_subdata<scalar_type>(world, boxes[me]);
    auto reference_data = get_subdata<scalar_type>(world, rotate_boxes[me]);
    auto output_data    = std::vector<scalar_type>(rotate_boxes[me].count());

    if (std::is_same<scalar_type, float>::value){
        // sometimes, run two tests to make sure there is no internal corruption
        // there is no need to do that for every data type
        reshape->apply(1, input_data.data(), output_data.data(), workspace.data());
        output_data = std::vector<scalar_type>(rotate_boxes[me].count());
        reshape->apply(1, input_data.data(), output_data.data(), workspace.data());
    }else{
        reshape->apply(1, input_data.data(), output_data.data(), workspace.data());
    }

//     mpi::dump(0, input_data,     "input");
//     mpi::dump(0, output_data,    "output");
//     mpi::dump(0, reference_data, "reference");

    tassert(match(output_data, reference_data));

    // do the same with a batch of 3
    int const batch_size = 4;
    std::vector<scalar_type> batch_input = input_data;
    std::vector<scalar_type> batch_reference = reference_data;
    for(int i=0; i<batch_size-1; i++){
        batch_input.insert(batch_input.end(), input_data.begin(), input_data.end());
        batch_reference.insert(batch_reference.end(), reference_data.begin(), reference_data.end());
        for(size_t j=batch_input.size() - input_data.size(); j < batch_input.size(); j++)
            batch_input[j] *= (i + 2);
        for(size_t j=batch_reference.size() - reference_data.size(); j < batch_reference.size(); j++)
            batch_reference[j] *= (i + 2);
    }
    output_data = std::vector<scalar_type>(batch_size * rotate_boxes[me].count());

    workspace.resize(batch_size * reshape->size_workspace());
    reshape->apply(batch_size, batch_input.data(), output_data.data(), workspace.data());

    tassert(match(output_data, batch_reference));
}

#ifdef Heffte_ENABLE_GPU
// splits the world box into a set of boxes with gird given by proc_grid
template<int hfast, int hmid, int hslow, int pfast, int pmid, int pslow, typename scalar_type, typename location_tag, reshape_algorithm variant>
void test_gpu(MPI_Comm const comm){
    /*
     * similar to the CPU case, but the data is located on the GPU
     */
    current_test<scalar_type, using_mpi, location_tag> test("-np " + std::to_string(mpi::comm_size(comm)) + "  "
                                                           + get_description<variant>(), comm);
    tassert( pfast * pmid * pslow == heffte::mpi::comm_size(comm) );

    int const me = heffte::mpi::comm_rank(comm);
    int const shift = 3;

    box3d<> world = {{0, 0, 0}, {hfast, hmid, hslow}};

    auto boxes   = split_world(world, {pfast, pmid, pslow});
    auto pencils = split_world(world, {pfast,    1, pmid * pslow});

    std::vector<box3d<>> rotate_boxes;
    if (std::is_same<scalar_type, std::complex<float>>::value){
        // shuffle the pencil boxes in some tests to check the case when there is no overlap between inbox and outbox
        // for the 2 by 2 grid, this shuffle ensures no overlap
        for(size_t i=0; i<boxes.size(); i++) rotate_boxes.push_back( pencils[(i + shift) % boxes.size()] );
    }else{
        for(auto b : pencils) rotate_boxes.push_back(b);
    }

    // create caches for a reshape algorithm, including creating a new mpi comm
    backend::device_instance<location_tag> device;
    auto reshape = make_test_reshape3d<location_tag, variant>(device.stream(), boxes, rotate_boxes, comm);
    gpu::vector<scalar_type> workspace(reshape->size_workspace());

    auto input_data     = get_subdata<scalar_type>(world, boxes[me]);
    auto cuinput_data   = gpu::transfer().load(input_data);
    auto reference_data = get_subdata<scalar_type>(world, rotate_boxes[me]);
    auto output_data    = gpu::vector<scalar_type>(rotate_boxes[me].count());

    if (std::is_same<scalar_type, float>::value){
        // sometimes, run two tests to make sure there is no internal corruption
        // there is no need to do that for every data type
        reshape->apply(1, cuinput_data.data(), output_data.data(), workspace.data());
        output_data = gpu::vector<scalar_type>(rotate_boxes[me].count());
        reshape->apply(1, cuinput_data.data(), output_data.data(), workspace.data());
    }else{
        reshape->apply(1, cuinput_data.data(), output_data.data(), workspace.data());
    }

    //auto ramout = cuda::unload(output_data);
    //mpi::dump(3, reference_data,     "reference_data");
    //mpi::dump(3, ramout,             "ramout");

    tassert(match(output_data, reference_data));
}
#endif

template<typename scalar_type>
void test_direct_reordered(MPI_Comm const comm){
    assert(mpi::comm_size(comm) == 4); // the rest is designed for 4 ranks
    current_test<scalar_type> test("-np " + std::to_string(mpi::comm_size(comm)) + "  direct_packer (unordered)", comm);

    box3d<> ordered_world(std::array<int, 3>{0, 0, 0}, std::array<int, 3>{8, 9, 10});

    auto ordered_inboxes  = split_world(ordered_world, {1, 2, 2});
    auto ordered_outboxes = split_world(ordered_world, {2, 2, 1});

    int const me = heffte::mpi::comm_rank(comm);
    auto input     = get_subdata<scalar_type>(ordered_world, ordered_inboxes[me]);
    auto reference = get_subdata<scalar_type>(ordered_world, ordered_outboxes[me]);

    box3d<> world(std::array<int, 3>{0, 0, 0}, std::array<int, 3>{9, 10, 8}, std::array<int, 3>{2, 0, 1});

    std::vector<box3d<>> inboxes, outboxes;
    for(auto b : split_world(world, {2, 2, 1})) inboxes.push_back({b.low, b.high, world.order});
    std::vector<box3d<>> temp = split_world(world, {2, 1, 2}); // need to swap the middle two entries
    for(auto i : std::vector<int>{0, 2, 1, 3}) outboxes.push_back({temp[i].low, temp[i].high, world.order});

    {
        auto reshape = make_reshape3d_alltoallv<tag::cpu>(nullptr, inboxes, outboxes, false, comm);
        std::vector<scalar_type> result(ordered_outboxes[me].count());
        std::vector<scalar_type> workspace(reshape->size_workspace());
        reshape->apply(1, input.data(), result.data(), workspace.data());

        tassert(match(result, reference));
    }{
        auto reshape = make_reshape3d_pointtopoint<tag::cpu>(nullptr, inboxes, outboxes, reshape_algorithm::p2p, false, comm);
        std::vector<scalar_type> result(ordered_outboxes[me].count());
        std::vector<scalar_type> workspace(reshape->size_workspace());
        reshape->apply(1, input.data(), result.data(), workspace.data());

        tassert(match(result, reference));
    }

    #ifdef Heffte_ENABLE_GPU
    for(bool use_gpu_aware : std::array<bool, 2>{true, false}){
    {
        backend::device_instance<tag::gpu> device;
        auto reshape = make_reshape3d_alltoallv<tag::gpu>(device.stream(), inboxes, outboxes, use_gpu_aware, comm);
        gpu::vector<scalar_type> workspace(reshape->size_workspace());

        auto cuinput = gpu::transfer().load(input);
        gpu::vector<scalar_type> curesult(ordered_outboxes[me].count());

        reshape->apply(1, cuinput.data(), curesult.data(), workspace.data());

        tassert(match(curesult, reference));
    }{
        backend::device_instance<tag::gpu> device;
        auto reshape = make_reshape3d_pointtopoint<tag::gpu>(device.stream(), inboxes, outboxes, reshape_algorithm::p2p_plined, use_gpu_aware, comm);
        gpu::vector<scalar_type> workspace(reshape->size_workspace());

        auto cuinput = gpu::transfer().load(input);
        gpu::vector<scalar_type> curesult(ordered_outboxes[me].count());

        reshape->apply(1, cuinput.data(), curesult.data(), workspace.data());
        tassert(match(curesult, reference));
    }}
    #endif
}

template<typename scalar_type, reshape_algorithm variant>
void test_reshape_transposed(MPI_Comm comm){
    assert(mpi::comm_size(comm) == 4); // the rest is designed for 4 ranks
    current_test<scalar_type> test("-np " + std::to_string(mpi::comm_size(comm)) + "  " + get_description<variant>() + " transposed", comm);

    box3d<> world(std::array<int, 3>{0, 0, 0}, std::array<int, 3>{1, 2, 3});

    std::vector<box3d<>> inboxes = split_world(world, {2, 1, 2});
    std::vector<box3d<>> ordered_outboxes = split_world(world, {1, 1, 4});

    int const me = heffte::mpi::comm_rank(comm);
    std::vector<scalar_type> input = get_subdata<scalar_type>(world, inboxes[me]);


    // the idea of the test is to try all combinations of order
    // test that MPI reshape with direct packer + on-node transpose reshape is equivalent to MPI reshape with transpose packer
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            if (i != j){
                int k = -1;
                for(int kk=0; kk<3; kk++) if (kk != i and kk != j) k = kk;
                std::array<int, 3> order = {i, j, k};
                if (i == 0 and j == 1 and k == 2) continue; // no transpose, no need to check

                std::vector<box3d<>> outboxes;
                for(auto b : ordered_outboxes) outboxes.push_back(box3d<>(b.low, b.high, order));

                heffte::plan_options options = default_options<default_cpu_backend>();
                options.algorithm = variant;

                auto mpi_tanspose_shaper  = make_reshape3d<default_cpu_backend>(nullptr, inboxes, outboxes, comm, options);
                auto mpi_direct_shaper    = make_reshape3d<default_cpu_backend>(nullptr, inboxes, ordered_outboxes, comm, options);
                auto cpu_transpose_shaper = make_reshape3d<default_cpu_backend>(nullptr, ordered_outboxes, outboxes, comm, options);

                std::vector<scalar_type> result(outboxes[me].count());
                std::vector<scalar_type> reference(outboxes[me].count());
                std::vector<scalar_type> workspace( // allocate one workspace vector for all reshape operations
                    std::max(std::max(mpi_tanspose_shaper->size_workspace(), mpi_direct_shaper->size_workspace()), cpu_transpose_shaper->size_workspace())
                );

                mpi_direct_shaper->apply(1, input.data(), result.data(), workspace.data());
                cpu_transpose_shaper->apply(1, result.data(), reference.data(), workspace.data());

                std::fill(result.begin(), result.end(), 0.0); // clear the temporary
                mpi_tanspose_shaper->apply(1, input.data(), result.data(), workspace.data());

                tassert(match(result, reference));

                #ifdef Heffte_ENABLE_GPU
                backend::device_instance<tag::gpu> device;
                heffte::plan_options cuoptions = default_options<gpu_backend>();
                cuoptions.algorithm = variant;

                auto cumpi_tanspose_shaper = make_reshape3d<gpu_backend>(device.stream(), inboxes, outboxes, comm, cuoptions);
                auto cumpi_direct_shaper   = make_reshape3d<gpu_backend>(device.stream(), inboxes, ordered_outboxes, comm, cuoptions);
                auto cuda_transpose_shaper = make_reshape3d<gpu_backend>(device.stream(), ordered_outboxes, outboxes, comm, cuoptions);

                gpu::vector<scalar_type> cuinput = gpu::transfer().load(input);
                gpu::vector<scalar_type> curesult(outboxes[me].count());
                gpu::vector<scalar_type> cureference(outboxes[me].count());
                gpu::vector<scalar_type> cuworkspace( // allocate one workspace vector for all reshape operations
                    std::max(std::max(cumpi_tanspose_shaper->size_workspace(), cumpi_direct_shaper->size_workspace()), cuda_transpose_shaper->size_workspace())
                );

                cumpi_direct_shaper->apply(1, cuinput.data(), curesult.data(), cuworkspace.data());
                cuda_transpose_shaper->apply(1, curesult.data(), cureference.data(), cuworkspace.data());

                curesult = gpu::vector<scalar_type>(outboxes[me].count());
                cumpi_tanspose_shaper->apply(1, cuinput.data(), curesult.data(), cuworkspace.data());

                tassert(match(curesult, gpu::transfer().unload(cureference)));
                #endif
            }
        }
    }
}

template<int hfast, int hmid, int hslow, int pfast, int pmid, int pslow, typename scalar_type, typename location_tag, reshape_algorithm variant>
void test_alltoone(MPI_Comm const comm){

    using backend_tag = typename default_tag<location_tag>::backend_tag;

    int const nprocs = mpi::comm_size(comm);
    int const me = mpi::comm_rank(comm);

    box3d<> const world = {{0, 0, 0}, {hfast, hmid, hslow}};
    box3d<> const empty_box = {{0, 0, 0}, {-1, -1, -1}};

    std::array<int,3> proc_i = {pfast, pmid, pslow};
    std::vector<box3d<>> inboxes  = heffte::split_world(world, proc_i);

    std::vector<box3d<>> outboxes;
    outboxes.push_back(world);
    for(int i=0; i<nprocs-1; i++) outboxes.push_back(empty_box);

    heffte::plan_options options = heffte::default_options<backend::stock>();
    options.algorithm = variant;
    options.use_reorder = false;

    backend::device_instance<location_tag> device;
    auto reshape = make_reshape3d<backend_tag>(device.stream(), inboxes, outboxes, comm, options);

    std::vector<scalar_type> world_vec(world.count());
    std::iota(world_vec.begin(), world_vec.end(), 1.0);

    std::vector<scalar_type> invec = get_subbox(world, inboxes[me], world_vec);

    auto loaded_invec = test_traits<location_tag>::load(invec);
    auto outvec = make_buffer_container<scalar_type>(device.stream(), reshape->size_output());
    auto workspace = make_buffer_container<scalar_type>(device.stream(), reshape->size_workspace());

    reshape->apply(1, loaded_invec.data(), outvec.data(), workspace.data());

    if (me == 0){
        tassert(approx(world_vec, outvec));
    }else{
        tassert(outvec.size() == 0);
    }
}

template<typename location_tag, reshape_algorithm variant>
void test_alltoone_variants(){
    MPI_Comm const comm = MPI_COMM_WORLD;

    current_test<float> test("-np " + std::to_string(mpi::comm_size(comm)) + "  " + get_description<variant>() + " all-2-1", comm);

    switch(mpi::comm_size(comm)) {
        case 4:
            test_alltoone<15, 15, 10, 2, 2, 1, float, location_tag, variant>(comm);
            test_alltoone<15, 15, 10, 1, 2, 2, float, location_tag, variant>(comm);
            test_alltoone<15, 15, 10, 1, 1, 4, float, location_tag, variant>(comm);
            test_alltoone<15, 1, 10, 1, 1, 4, float, location_tag, variant>(comm);
            test_alltoone<15, 15, 10, 1, 1, 4, float, location_tag, variant>(comm);
            test_alltoone<1, 15, 10, 1, 1, 4, float, location_tag, variant>(comm);
            test_alltoone<1, 15, 10, 1, 2, 2, float, location_tag, variant>(comm);
            test_alltoone<15, 15, 10, 2, 2, 1, float, location_tag, variant>(comm);
            break;
        case 12:
            test_alltoone<15, 15, 10, 3, 4, 1, float, location_tag, variant>(comm);
            test_alltoone<15, 25, 10, 1, 6, 2, float, location_tag, variant>(comm);
            test_alltoone<25, 25, 10, 1, 3, 4, float, location_tag, variant>(comm);
            test_alltoone<15, 1, 10, 3, 1, 4, float, location_tag, variant>(comm);
            test_alltoone<15, 25, 10, 4, 1, 3, float, location_tag, variant>(comm);
            test_alltoone<1, 15, 21, 1, 1, 12, float, location_tag, variant>(comm);
            test_alltoone<1, 15, 10, 1, 3, 4, float, location_tag, variant>(comm);
            test_alltoone<25, 15, 10, 4, 3, 1, float, location_tag, variant>(comm);
            test_alltoone<9, 11, 5, 2, 3, 2, float, location_tag, variant>(comm);
            test_alltoone<9, 11, 5, 1, 4, 3, float, location_tag, variant>(comm);
        default:
            break;
    }
}

void test_alltoone_all(){
    test_alltoone_variants<tag::cpu, reshape_algorithm::alltoall>();
    test_alltoone_variants<tag::cpu, reshape_algorithm::alltoallv>();
    test_alltoone_variants<tag::cpu, reshape_algorithm::p2p>();
    test_alltoone_variants<tag::cpu, reshape_algorithm::p2p_plined>();
    #ifdef Heffte_ENABLE_GPU
    test_alltoone_variants<tag::gpu, reshape_algorithm::alltoall>();
    test_alltoone_variants<tag::gpu, reshape_algorithm::alltoallv>();
    test_alltoone_variants<tag::gpu, reshape_algorithm::p2p>();
    test_alltoone_variants<tag::gpu, reshape_algorithm::p2p_plined>();
    #endif
}

void perform_tests_cpu(){
    MPI_Comm const comm = MPI_COMM_WORLD;

    test_boxes(comm);

    switch(mpi::comm_size(comm)) {
        // note that the number of boxes must match the comm size
        // that is the product of the last three of the box dimensions
        case 4:
            test_cpu<10, 13, 10, 2, 2, 1, float, heffte::tag::cpu, reshape_algorithm::alltoall>(comm);
            test_cpu<10, 20, 17, 2, 2, 1, double, heffte::tag::cpu, reshape_algorithm::alltoall>(comm);
            test_cpu<30, 10, 10, 2, 2, 1, std::complex<float>, heffte::tag::cpu, reshape_algorithm::alltoall>(comm);
            test_cpu<11, 10, 13, 2, 2, 1, std::complex<double>, heffte::tag::cpu, reshape_algorithm::alltoall>(comm);
            test_cpu<10, 13, 10, 2, 2, 1, float, heffte::tag::cpu, reshape_algorithm::alltoallv>(comm);
            test_cpu<10, 20, 17, 2, 2, 1, double, heffte::tag::cpu, reshape_algorithm::alltoallv>(comm);
            test_cpu<30, 10, 10, 2, 2, 1, std::complex<float>, heffte::tag::cpu, reshape_algorithm::alltoallv>(comm);
            test_cpu<11, 10, 13, 2, 2, 1, std::complex<double>, heffte::tag::cpu, reshape_algorithm::alltoallv>(comm);
            test_cpu<10, 13, 10, 2, 2, 1, float, heffte::tag::cpu, reshape_algorithm::p2p>(comm);
            test_cpu<10, 20, 17, 2, 2, 1, double, heffte::tag::cpu, reshape_algorithm::p2p>(comm);
            test_cpu<30, 10, 10, 2, 2, 1, std::complex<float>, heffte::tag::cpu, reshape_algorithm::p2p>(comm);
            test_cpu<11, 10, 13, 2, 2, 1, std::complex<double>, heffte::tag::cpu, reshape_algorithm::p2p>(comm);
            break;
        case 12:
            test_cpu<13, 13, 10, 3, 4, 1, float, heffte::tag::cpu, reshape_algorithm::alltoall>(comm);
            test_cpu<41, 17, 15, 3, 2, 2, std::complex<double>, heffte::tag::cpu, reshape_algorithm::alltoall>(comm);
            test_cpu<13, 13, 10, 3, 4, 1, float, heffte::tag::cpu, reshape_algorithm::alltoallv>(comm);
            test_cpu<16, 21, 17, 2, 3, 2, double, heffte::tag::cpu, reshape_algorithm::alltoallv>(comm);
            test_cpu<38, 13, 20, 1, 4, 3, std::complex<float>, heffte::tag::cpu, reshape_algorithm::alltoallv>(comm);
            test_cpu<41, 17, 15, 3, 2, 2, std::complex<double>, heffte::tag::cpu, reshape_algorithm::alltoallv>(comm);
            test_cpu<13, 13, 10, 3, 4, 1, float, heffte::tag::cpu, reshape_algorithm::p2p>(comm);
            test_cpu<16, 21, 17, 2, 3, 2, double, heffte::tag::cpu, reshape_algorithm::p2p>(comm);
            test_cpu<38, 13, 20, 1, 4, 3, std::complex<float>, heffte::tag::cpu, reshape_algorithm::p2p>(comm);
            test_cpu<41, 17, 15, 3, 2, 2, std::complex<double>, heffte::tag::cpu, reshape_algorithm::p2p>(comm);
            break;
        default:
            // unknown test
            break;
    }
}

void perform_tests_gpu(){
    #ifdef Heffte_ENABLE_GPU
    MPI_Comm const comm = MPI_COMM_WORLD;

    switch(mpi::comm_size(comm)) {
        // note that the number of boxes must match the comm size
        // that is the product of the last three of the box dimensions
        case 4:
            test_gpu<10, 13, 10, 2, 2, 1, float, heffte::tag::gpu, reshape_algorithm::alltoall>(comm);
            test_gpu<10, 20, 17, 2, 2, 1, double, heffte::tag::gpu, reshape_algorithm::alltoall>(comm);
            test_gpu<30, 10, 10, 2, 2, 1, std::complex<float>, heffte::tag::gpu, reshape_algorithm::alltoall>(comm);
            test_gpu<11, 10, 13, 2, 2, 1, std::complex<double>, heffte::tag::gpu, reshape_algorithm::alltoall>(comm);
            test_gpu<10, 13, 10, 2, 2, 1, float, heffte::tag::gpu, reshape_algorithm::alltoallv>(comm);
            test_gpu<10, 20, 17, 2, 2, 1, double, heffte::tag::gpu, reshape_algorithm::alltoallv>(comm);
            test_gpu<30, 10, 10, 2, 2, 1, std::complex<float>, heffte::tag::gpu, reshape_algorithm::alltoallv>(comm);
            test_gpu<11, 10, 13, 2, 2, 1, std::complex<double>, heffte::tag::gpu, reshape_algorithm::alltoallv>(comm);
            test_gpu<10, 13, 10, 2, 2, 1, float, heffte::tag::gpu, reshape_algorithm::p2p>(comm);
            test_gpu<10, 20, 17, 2, 2, 1, double, heffte::tag::gpu, reshape_algorithm::p2p>(comm);
            test_gpu<30, 10, 10, 2, 2, 1, std::complex<float>, heffte::tag::gpu, reshape_algorithm::p2p>(comm);
            test_gpu<11, 10, 13, 2, 2, 1, std::complex<double>, heffte::tag::gpu, reshape_algorithm::p2p>(comm);
            break;
        case 12:
            test_gpu<13, 13, 10, 3, 4, 1, float, heffte::tag::gpu, reshape_algorithm::alltoall>(comm);
            test_gpu<41, 17, 15, 3, 2, 2, std::complex<double>, heffte::tag::gpu, reshape_algorithm::alltoall>(comm);
            test_gpu<13, 13, 10, 3, 4, 1, float, heffte::tag::gpu, reshape_algorithm::alltoallv>(comm);
            test_gpu<16, 21, 17, 2, 3, 2, double, heffte::tag::gpu, reshape_algorithm::alltoallv>(comm);
            test_gpu<38, 13, 20, 1, 4, 3, std::complex<float>, heffte::tag::gpu, reshape_algorithm::alltoallv>(comm);
            test_gpu<41, 17, 15, 3, 2, 2, std::complex<double>, heffte::tag::gpu, reshape_algorithm::alltoallv>(comm);
            test_gpu<13, 13, 10, 3, 4, 1, float, heffte::tag::gpu, reshape_algorithm::p2p>(comm);
            test_gpu<16, 21, 17, 2, 3, 2, double, heffte::tag::gpu, reshape_algorithm::p2p>(comm);
            test_gpu<38, 13, 20, 1, 4, 3, std::complex<float>, heffte::tag::gpu, reshape_algorithm::p2p>(comm);
            test_gpu<41, 17, 15, 3, 2, 2, std::complex<double>, heffte::tag::gpu, reshape_algorithm::p2p>(comm);
            break;
        default:
            // unknown test
            break;
    }
    #endif
}

void perform_tests_reorder(){
    MPI_Comm const comm = MPI_COMM_WORLD;

    if (mpi::comm_size(comm) == 4){
        test_direct_reordered<double>(comm);
        test_direct_reordered<std::complex<float>>(comm);
        test_reshape_transposed<float, reshape_algorithm::alltoallv>(comm);
        test_reshape_transposed<std::complex<double>, reshape_algorithm::alltoallv>(comm);
        test_reshape_transposed<float, reshape_algorithm::alltoall>(comm);
        test_reshape_transposed<std::complex<double>, reshape_algorithm::alltoall>(comm);
        test_reshape_transposed<double, reshape_algorithm::p2p>(comm);
        test_reshape_transposed<std::complex<float>, reshape_algorithm::p2p>(comm);
        test_reshape_transposed<double, reshape_algorithm::p2p_plined>(comm);
        test_reshape_transposed<std::complex<float>, reshape_algorithm::p2p_plined>(comm);
    }
}

void perform_all_tests(){
    all_tests<> name("heffte reshape methods");
    perform_tests_cpu();
    perform_tests_gpu();
    perform_tests_reorder();
    test_alltoone_all();
}

int main(int argc, char *argv[]){

    MPI_Init(&argc, &argv);

    perform_all_tests();

    MPI_Finalize();

    return 0;
}
