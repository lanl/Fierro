/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
       Performance test for 3D FFTs using heFFTe
*/

#include "test_fft3d.h"

#ifdef BENCH_C2C
template<typename precision>
struct bench_types{
    using input = std::complex<precision>;
    using output = std::complex<precision>;
};
#endif
#ifdef BENCH_R2C
template<typename precision>
struct bench_types{
    using input = precision;
    using output = std::complex<precision>;
};
#endif
#ifdef BENCH_R2R
template<typename precision>
struct bench_types{
    using input = precision;
    using output = precision;
};
#endif

template<typename backend_tag, typename precision_type, typename index>
void benchmark_fft(std::array<int,3> size_fft, std::deque<std::string> const &args){

    using input_type = typename bench_types<precision_type>::input;
    using output_type = typename bench_types<precision_type>::output;

    int me, nprocs;
    MPI_Comm fft_comm = MPI_COMM_WORLD;  // Change if need to compute FFT within a subcommunicator
    MPI_Comm_rank(fft_comm, &me);
    MPI_Comm_size(fft_comm, &nprocs);

    // Create input and output boxes on local processor
    box3d<index> const world = {{0, 0, 0}, {size_fft[0]-1, size_fft[1]-1, size_fft[2]-1}};

    // Get grid of processors at input and output
    std::array<int,3> proc_i = heffte::proc_setup_min_surface(world, nprocs);
    std::array<int,3> proc_o = heffte::proc_setup_min_surface(world, nprocs);

    // Check if user in/out processor grids are pencil-shaped, useful for performance comparison with other libraries
    if (has_option(args, "-io_pencils")){
        std::array<int, 2> proc_grid = make_procgrid(nprocs);
        proc_i = {1, proc_grid[0], proc_grid[1]};
        proc_o = {1, proc_grid[0], proc_grid[1]};
    }

    if (has_option(args, "-ingrid"))
        proc_i = get_grid(args, "-ingrid");
    if (has_option(args, "-outgrid"))
        proc_o = get_grid(args, "-outgrid");

    if (proc_i[0] * proc_i[1] * proc_i[2] != nprocs)
        throw std::runtime_error(std::string(" incorrect input processor grid, specified ") + std::to_string(proc_i[0] * proc_i[1] * proc_i[2]) + " processors but using " + std::to_string(nprocs));
    if (proc_o[0] * proc_o[1] * proc_o[2] != nprocs)
        throw std::runtime_error(std::string(" incorrect output processor grid, specified ") + std::to_string(proc_o[0] * proc_o[1] * proc_o[2]) + " processors but using " + std::to_string(nprocs));

    int const r2c_dir = [&]()->int{
        #ifdef BENCH_R2C
        int r = get_int_arg("-r2c_dir", args);
        return (r == -1) ? 0 : r;
        #else
        return -1;
        #endif
    }();

    int const batch_size = get_int_arg("-batch", args, 1);

    std::vector<box3d<index>> inboxes  = heffte::split_world(world, proc_i);
    std::vector<box3d<index>> outboxes = heffte::split_world((r2c_dir == -1) ? world : world.r2c(r2c_dir), proc_o);

    #ifdef Heffte_ENABLE_GPU
    if (backend::uses_gpu<backend_tag>::value and has_option(args, "-mps")){
        heffte::gpu::device_set(me % heffte::gpu::device_count());
    }
    #endif

    // Define 3D FFT plan

    heffte::plan_options options = args_to_options<backend_tag>(args);

    // the call above uses the following plan, get it twice to give verbose info of the grid-shapes
    logic_plan3d<index> plan = plan_operations<index>({inboxes, outboxes}, r2c_dir, options, me);

    // computes the grids
    std::vector<std::array<int, 3>> grids = compute_grids(plan);

    // writes out the proc_grid in the given iteration, skips if no reshape was applied
    auto print_proc_grid = [&](int i){
        if (i > 0 and grids[i] == grids[i-1]){  cout << "-  ";  return;  }
        cout << "(" << grids[i][0] << ", " << grids[i][1] << ", " << grids[i][2] << ")  ";
    };

    // Locally initialize input
    auto input = make_data<input_type>(batch_size, inboxes[me]);
    auto reference_input = input; // safe a copy for error checking

    // define allocation for in-place transform
    #if defined(BENCH_C2C) || defined(BENCH_R2R)
    auto fft = make_fft3d<backend_tag>(inboxes[me], outboxes[me], fft_comm, options);
    std::vector<output_type> output(batch_size * std::max(fft.size_outbox(), fft.size_inbox()));
    std::copy(input.begin(), input.end(), output.begin());

    output_type *output_array = output.data();
    #ifdef Heffte_ENABLE_GPU
    gpu::vector<output_type> gpu_output;
    if (backend::uses_gpu<backend_tag>::value){
        gpu_output = gpu::transfer::load(output);
        output_array = gpu_output.data();
    }
    #endif
    output_type *input_array = output_array; // using in-place transform for the C2C case
    output_type *result = output.data();
    #endif

    #ifdef BENCH_R2C
    auto fft = make_fft3d_r2c<backend_tag>(inboxes[me], outboxes[me], r2c_dir, fft_comm, options);

    input_type *input_array = input.data();
    #ifdef Heffte_ENABLE_GPU
    gpu::vector<input_type> gpu_input;
    if (backend::uses_gpu<backend_tag>::value){
        gpu_input = gpu::transfer::load(input);
        input_array = gpu_input.data();
    }
    #endif

    typename heffte::fft3d<backend_tag>::template buffer_container<output_type> output(batch_size * fft.size_outbox());
    output_type *output_array = output.data();
    input_type *result = input.data();
    #endif

    // Define workspace array
    typename heffte::fft3d<backend_tag>::template buffer_container<output_type> workspace(batch_size * fft.size_workspace());

    // Warmup
    heffte::add_trace("mark warmup begin");
    fft.forward(input_array, output_array,  scale::full);
    fft.backward(output_array, input_array);

    // Execution
    int const ntest = nruns(args);
    MPI_Barrier(fft_comm);
    double t = -MPI_Wtime();
    if (batch_size == 1){
        for(int i = 0; i < ntest; ++i) {
            heffte::add_trace("mark forward begin");
            fft.forward(input_array, output_array, workspace.data(), scale::full);
            heffte::add_trace("mark backward begin");
            fft.backward(output_array, input_array, workspace.data());
        }
    }else{
        for(int i = 0; i < ntest; ++i) {
            heffte::add_trace("mark forward begin with batch: " + std::to_string(batch_size));
            fft.forward(batch_size, input_array, output_array, workspace.data(), scale::full);
            heffte::add_trace("mark backward begin");
            fft.backward(batch_size, output_array, input_array, workspace.data());
        }
    }
    #ifdef Heffte_ENABLE_GPU
    if (backend::uses_gpu<backend_tag>::value)
        gpu::synchronize_default_stream();
    #endif
    MPI_Barrier(fft_comm);
    t += MPI_Wtime();

    // Get execution time
    double t_max = 0.0;
	MPI_Reduce(&t, &t_max, 1, MPI_DOUBLE, MPI_MAX, 0, fft_comm);

    // Validate result, but first unload from the GPU
    #ifdef Heffte_ENABLE_GPU
    #if defined(BENCH_C2C) || defined(BENCH_R2R)
    if (std::is_same<backend_tag, gpu_backend>::value){
        output = gpu::transfer::unload(gpu_output);
        result = output.data();
    }
    #endif
    #ifdef BENCH_R2C
    if (std::is_same<backend_tag, gpu_backend>::value){
        input = gpu::transfer::unload(gpu_input);
        result = input.data();
    }
    #endif
    #endif

    precision_type err = 0.0;
    for(size_t i=0; i<input.size(); i++)
        err = std::max(err, std::abs(reference_input[i] - result[i]));
    precision_type mpi_max_err = 0.0;
    MPI_Allreduce(&err, &mpi_max_err, 1, mpi::type_from<precision_type>(), MPI_MAX, fft_comm);

    if (mpi_max_err > precision<std::complex<precision_type>>::tolerance){
        // benchmark failed, the error is too much
        if (me == 0){
            cout << "------------------------------- \n"
                 << "ERROR: observed error after heFFTe benchmark exceeds the tolerance\n"
                 << "       tolerance: " << precision<std::complex<precision_type>>::tolerance
                 << "  error: " << mpi_max_err << endl;
        }
        return;
    }

    // Print results
    if(me==0){
        t_max = t_max / (2.0 * ntest);
        double const fftsize  = static_cast<double>(world.count());
        double const floprate = 5.0 * batch_size * fftsize * std::log(fftsize) * 1e-9 / std::log(2.0) / t_max;
        long long mem_usage = static_cast<long long>(fft.size_inbox()) + static_cast<long long>(fft.size_outbox())
                            + static_cast<long long>(fft.size_workspace());
        mem_usage *= sizeof(output_type);
        mem_usage /= 1024ll * 1024ll; // convert to MB
        cout << "\n----------------------------------------------------------------------------- \n";
        cout << "heFFTe performance test\n";
        cout << "----------------------------------------------------------------------------- \n";
        cout << "Backend:   " << backend::name<backend_tag>() << "\n";
        cout << "Size:      " << world.size[0] << "x" << world.size[1] << "x" << world.size[2] << "\n";
        cout << "MPI ranks: " << setw(4) << nprocs << "\n";
        cout << "Grids: ";
        for(int i=0; i<5; i++)
            print_proc_grid(i);
        cout << "\n";
        cout << "Time per run: " << t_max << " (s)\n";
        cout << "Performance:  " << floprate << " GFlops/s\n";
        cout << "Memory usage: " << mem_usage << "MB/rank\n";
        cout << "Tolerance:    " << precision<std::complex<precision_type>>::tolerance << "\n";
        cout << "Max error:    " << mpi_max_err << "\n";
        cout << endl;
    }
}

template<typename backend_tag>
bool perform_benchmark(std::string const &precision_string, std::string const &backend_string, std::string const &backend_name,
                       std::array<int,3> size_fft, std::deque<std::string> const &args){
    if (backend_string == backend_name){
        if (precision_string == "float"){
            benchmark_fft<backend_tag, float, int>(size_fft, args);
        }else if (precision_string == "double"){
            benchmark_fft<backend_tag, double, int>(size_fft, args);
        }else if (precision_string == "float-long"){
            benchmark_fft<backend_tag, float, long long>(size_fft, args);
        }else{ // double-long
            benchmark_fft<backend_tag, double, long long>(size_fft, args);
        }
        return true;
    }
    return false;
}

int main(int argc, char *argv[]){

    MPI_Init(&argc, &argv);

    #ifdef BENCH_C2C
    std::string bench_executable = "./speed3d_c2c";
    #endif
    #ifdef BENCH_R2C
    std::string bench_executable = "./speed3d_r2c";
    #endif
    #ifdef BENCH_R2R
    std::string bench_executable = "./speed3d_r2r";
    #endif

#ifdef BENCH_R2R
    std::string backends = "stock-cos stock-sin ";
    #ifdef Heffte_ENABLE_FFTW
    backends += "fftw-cos fftw-sin ";
    #endif
    #ifdef Heffte_ENABLE_CUDA
    backends += "cufft-cos cufft-sin ";
    #endif
    #ifdef Heffte_ENABLE_ROCM
    backends += "rocfft-cos rocfft-sin ";
    #endif
    #ifdef Heffte_ENABLE_ONEAPI
    backends += "onemkl-cos onemkl-sin ";
    #endif
    #ifdef Heffte_ENABLE_MKL
    backends += "mkl-cos mkl-sin ";
    #endif
#else
    std::string backends = "stock ";
    #ifdef Heffte_ENABLE_FFTW
    backends += "fftw ";
    #endif
    #ifdef Heffte_ENABLE_CUDA
    backends += "cufft ";
    #endif
    #ifdef Heffte_ENABLE_ROCM
    backends += "rocfft ";
    #endif
    #ifdef Heffte_ENABLE_ONEAPI
    backends += "onemkl ";
    #endif
    #ifdef Heffte_ENABLE_MKL
    backends += "mkl ";
    #endif
#endif

    if (argc < 6){
        if (mpi::world_rank(0)){
            cout << "\nUsage:\n    mpirun -np x " << bench_executable << " <backend> <precision> <size-x> <size-y> <size-z> <args>\n\n"
                 << "    options\n"
                 << "        backend is the 1-D FFT library\n"
                 << "            available options for this build: " << backends << "\n"
                 << "        precision is either float or double\n"
                 << "          use float-long or double-long to enable 64-bit indexing\n"
                 << "        size-x/y/z are the 3D array dimensions \n\n"
                 << "        args is a set of optional arguments that define algorithmic tweaks and variations\n"
                 << "         -reorder: reorder the elements of the arrays so that each 1-D FFT will use contiguous data\n"
                 << "         -no-reorder: some of the 1-D will be strided (non contiguous)\n"
                 << "         -a2a: use MPI_Alltoall() communication method\n"
                 << "         -a2av: use MPI_Alltoallv() communication method (default)\n"
                 << "         -p2p: use MPI_Send() and MPI_Irecv() communication methods\n"
                 << "         -p2p_pl: use MPI_Isend() and MPI_Irecv() communication methods\n"
                 << "         -no-gpu-aware: move the data to the cpu before doing gpu operations (gpu backends only)\n"
                 << "         -pencils: use pencil reshape logic\n"
                 << "         -slabs: use slab reshape logic\n"
                 << "         -io_pencils: if input and output proc grids are pencils, useful for comparison with other libraries \n"
                 << "         -ingrid x y z: specifies the processor grid to use in the input, x y z must be integers \n"
                 << "         -outgrid x y z: specifies the processor grid to use in the output, x y z must be integers \n"
                 << "         -subcomm num_ranks: specifies the number of ranks to use in intermediate reshapes\n"
                 << "         -batch batch_size: specifies the size of the batch to use in the benchmark\n"
                 << "         -r2c_dir dir: specifies the r2c direction for the r2c tests, dir must be 0 1 or 2 \n"
                 << "         -mps: for the cufft backend and multiple gpus, associate the mpi ranks with different cuda devices\n"
                 << "         -nX: number of times to repeat the run, accepted variants are -n5 (default), -n10, -n50\n"
                 #ifdef BENCH_R2R
                 << "Examples:\n"
                 << "    mpirun -np  4 " << bench_executable << " fftw-cos  double 128 128 128 -p2p\n"
                 << "    mpirun -np  8 " << bench_executable << " cufft-cos float  256 256 256\n"
                 << "    mpirun -np 12 " << bench_executable << " fftw-sin  double 512 512 512 -slabs\n\n";
                 #else
                 << "Examples:\n"
                 << "    mpirun -np  4 " << bench_executable << " fftw  double 128 128 128 -no-reorder\n"
                 << "    mpirun -np  8 " << bench_executable << " cufft float  256 256 256\n"
                 << "    mpirun -np 12 " << bench_executable << " fftw  double 512 512 512 -p2p -slabs\n\n";
                 #endif
        }

        MPI_Finalize();
        return 0;
    }

    std::array<int,3> size_fft = { 0, 0, 0 };

    std::string backend_string = argv[1];

    std::string precision_string = argv[2];
    if (precision_string != "float"      and precision_string != "double" and
        precision_string != "float-long" and precision_string != "double-long"){
        if (mpi::world_rank(0)){
            std::cout << "Invalid precision!\n";
            std::cout << "Must use float or double" << std::endl;
        }
        MPI_Finalize();
        return 0;
    }

    try{
        size_fft = { std::stoi(argv[3]), std::stoi(argv[4]), std::stoi(argv[5])};
        for(auto s : size_fft) if (s < 1) throw std::invalid_argument("negative input");
    }catch(std::invalid_argument &e){
        if (mpi::world_rank(0)){
            std::cout << "Cannot convert the sizes into positive integers!\n";
            std::cout << "Encountered error: " << e.what() << std::endl;
        }
        MPI_Finalize();
        return 0;
    }

    init_tracing(bench_executable + "_" + backend_string + "_" + precision_string
                 + std::string(argv[3]) + "_" + std::string(argv[4]) + "_" + std::string(argv[5]));

    bool valid_backend = false;
#ifdef BENCH_R2R
    #ifdef Heffte_ENABLE_FFTW
    valid_backend = valid_backend or perform_benchmark<backend::fftw_cos>(precision_string, backend_string, "fftw-cos", size_fft, arguments(argc, argv));
    valid_backend = valid_backend or perform_benchmark<backend::fftw_sin>(precision_string, backend_string, "fftw-sin", size_fft, arguments(argc, argv));
    #endif
    valid_backend = valid_backend or perform_benchmark<backend::stock_cos>(precision_string, backend_string, "stock-cos", size_fft, arguments(argc, argv));
    valid_backend = valid_backend or perform_benchmark<backend::stock_sin>(precision_string, backend_string, "stock-sin", size_fft, arguments(argc, argv));
    #ifdef Heffte_ENABLE_MKL
    valid_backend = valid_backend or perform_benchmark<backend::mkl_cos>(precision_string, backend_string, "mkl-cos", size_fft, arguments(argc, argv));
    valid_backend = valid_backend or perform_benchmark<backend::mkl_sin>(precision_string, backend_string, "mkl-sin", size_fft, arguments(argc, argv));
    #endif
    #ifdef Heffte_ENABLE_CUDA
    valid_backend = valid_backend or perform_benchmark<backend::cufft_cos>(precision_string, backend_string, "cufft-cos", size_fft, arguments(argc, argv));
    valid_backend = valid_backend or perform_benchmark<backend::cufft_sin>(precision_string, backend_string, "cufft-sin", size_fft, arguments(argc, argv));
    #endif
    #ifdef Heffte_ENABLE_ROCM
    valid_backend = valid_backend or perform_benchmark<backend::rocfft_cos>(precision_string, backend_string, "rocfft-cos", size_fft, arguments(argc, argv));
    valid_backend = valid_backend or perform_benchmark<backend::rocfft_sin>(precision_string, backend_string, "rocfft-sin", size_fft, arguments(argc, argv));
    #endif
    #ifdef Heffte_ENABLE_ONEAPI
    valid_backend = valid_backend or perform_benchmark<backend::onemkl_cos>(precision_string, backend_string, "onemkl-cos", size_fft, arguments(argc, argv));
    valid_backend = valid_backend or perform_benchmark<backend::onemkl_sin>(precision_string, backend_string, "onemkl-sin", size_fft, arguments(argc, argv));
    #endif
#else
    #ifdef Heffte_ENABLE_FFTW
    valid_backend = valid_backend or perform_benchmark<backend::fftw>(precision_string, backend_string, "fftw", size_fft, arguments(argc, argv));
    #endif
    valid_backend = valid_backend or perform_benchmark<backend::stock>(precision_string, backend_string, "stock", size_fft, arguments(argc, argv));
    #ifdef Heffte_ENABLE_MKL
    valid_backend = valid_backend or perform_benchmark<backend::mkl>(precision_string, backend_string, "mkl", size_fft, arguments(argc, argv));
    #endif
    #ifdef Heffte_ENABLE_CUDA
    valid_backend = valid_backend or perform_benchmark<backend::cufft>(precision_string, backend_string, "cufft", size_fft, arguments(argc, argv));
    #endif
    #ifdef Heffte_ENABLE_ROCM
    valid_backend = valid_backend or perform_benchmark<backend::rocfft>(precision_string, backend_string, "rocfft", size_fft, arguments(argc, argv));
    #endif
    #ifdef Heffte_ENABLE_ONEAPI
    valid_backend = valid_backend or perform_benchmark<backend::onemkl>(precision_string, backend_string, "onemkl", size_fft, arguments(argc, argv));
    #endif
#endif

    if (not valid_backend){
        if (mpi::world_rank(0)){
            std::cout << "Invalid backend " << backend_string << "\n";
            std::cout << "The available backends are: " << backends << std::endl;
        }
        MPI_Finalize();
        return 0;
    }

    finalize_tracing();

    MPI_Finalize();
    return 0;
}
