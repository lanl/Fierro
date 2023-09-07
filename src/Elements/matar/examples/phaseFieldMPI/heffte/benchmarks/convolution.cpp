/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
       Simple convolution test for 3D FFTs using heFFTe
*/

// Calculates FFT^-1( FFT(X) .* FFT(Y) ), valid for the 3 CPU kernels
#include "test_fft3d.h"

#define BENCH_INPUT std::complex<precision_type>

template<typename backend_tag, typename precision_type, typename index>
void benchmark_convolution(std::array<int,3> size_fft, std::deque<std::string> const &args){

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

    std::vector<box3d<index>> inboxes  = heffte::split_world(world, proc_i);
    std::vector<box3d<index>> outboxes = heffte::split_world(world, proc_o);


    // Define 3D FFT plan
    heffte::plan_options options = args_to_options<backend_tag>(args);

    auto fft = make_fft3d<backend_tag>(inboxes[me], outboxes[me], fft_comm, options);

    std::array<int, 2> proc_grid = make_procgrid(nprocs);
    // writes out the proc_grid in the given dimension
    auto print_proc_grid = [&](int i){
        switch(i){
            case -1: cout << "(" << proc_i[0] << ", " << proc_i[1] << ", " << proc_i[2] << ")  "; break;
            case  0: cout << "(" << 1 << ", " << proc_grid[0] << ", " << proc_grid[1] << ")  "; break;
            case  1: cout << "(" << proc_grid[0] << ", " << 1 << ", " << proc_grid[1] << ")  "; break;
            case  2: cout << "(" << proc_grid[0] << ", " << proc_grid[1] << ", " << 1 << ")  "; break;
            case  3: cout << "(" << proc_o[0] << ", " << proc_o[1] << ", " << proc_o[2] << ")  "; break;
            default:
                throw std::runtime_error("printing incorrect direction");
        }
    };

    // the call above uses the following plan, get it twice to give verbose info of the grid-shapes
    logic_plan3d<index> plan = plan_operations<index>({inboxes, outboxes}, -1, heffte::default_options<backend_tag>(), me);

    // Locally initialize inputs
    auto X = make_data<BENCH_INPUT>(inboxes[me]);
    auto Y = X;

    // define allocation for in-place transform
    std::vector<std::complex<precision_type>> output(std::max(fft.size_outbox(), fft.size_inbox()));
    std::copy(X.begin(), X.end(), output.begin());

    std::complex<precision_type> *output_array = output.data();

    // Define workspace array
    typename heffte::fft3d<backend_tag>::template buffer_container<std::complex<precision_type>> workspace(fft.size_workspace());


    // Execution
    int const ntest = nruns(args);
    MPI_Barrier(fft_comm);
    double t = -MPI_Wtime();
        fft.forward(output_array, output_array, workspace.data(), scale::full);

        for(size_t i=0; i<output.size(); ++i){
            // cout << output_array[i] << endl;
            output_array[i] *= output_array[i];
        }
        fft.backward(output_array, output_array, workspace.data());

    MPI_Barrier(fft_comm);
    t += MPI_Wtime();

    // Get execution time
    double t_max = 0.0;
	MPI_Reduce(&t, &t_max, 1, MPI_DOUBLE, MPI_MAX, 0, fft_comm);

    // Print results
    if(me==0){
        t_max = t_max / (2.0 * ntest);
//         double const fftsize  = static_cast<double>(world.count());
//         double const floprate = 5.0 * fftsize * std::log(fftsize) * 1e-9 / std::log(2.0) / t_max;
//         long long mem_usage = static_cast<long long>(fft.size_inbox()) + static_cast<long long>(fft.size_outbox())
//                             + static_cast<long long>(fft.size_workspace());
//         mem_usage *= sizeof(std::complex<precision_type>);
//         mem_usage /= 1024ll * 1024ll; // convert to MB
        cout << "\n----------------------------------------------------------------------------- \n";
        cout << "heFFTe performance test\n";
        cout << "----------------------------------------------------------------------------- \n";
        cout << "Backend:   Convolution FFT^{-1} ( FFT(X).* FFT(Y) ) \n";
        cout << "Size input X:      " << world.size[0] << "x" << world.size[1] << "x" << world.size[2] << "\n";
        cout << "Size input Y:      " << world.size[0] << "x" << world.size[1] << "x" << world.size[2] << "\n";
        cout << "MPI ranks: " << setw(4) << nprocs << "\n";
        cout << "Grids: ";
        print_proc_grid(-1); // prints the initial grid
        for(int i=0; i<3; i++) // if reshape was applied, show the new intermediate grid assuming pencil reshape
            if (not match(plan.in_shape[i], plan.out_shape[i]) and not match(plan.out_shape[i], plan.out_shape[3])) print_proc_grid(i);
        print_proc_grid(3); // print the final grid
        cout << "\n";
        cout << "Time per run: " << t_max << " (s)\n";
        cout << endl;
    }
}

template<typename backend_tag>
bool perform_benchmark(std::string const &precision_string, std::string const &backend_string, std::string const &backend_name,
                       std::array<int,3> size_fft, std::deque<std::string> const &args){
    if (backend_string == backend_name){
        if (precision_string == "float"){
            benchmark_convolution<backend_tag, float, int>(size_fft, args);
        }else if (precision_string == "double"){
            benchmark_convolution<backend_tag, double, int>(size_fft, args);
        }else if (precision_string == "float-long"){
            benchmark_convolution<backend_tag, float, long long>(size_fft, args);
        }else{ // double-long
            benchmark_convolution<backend_tag, double, long long>(size_fft, args);
        }
        return true;
    }
    return false;
}

int main(int argc, char *argv[]){

    MPI_Init(&argc, &argv);

    std::string backends = "stock ";
    #ifdef Heffte_ENABLE_FFTW
    backends += "fftw ";
    #endif
    #ifdef Heffte_ENABLE_MKL
    backends += "mkl ";
    #endif

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

    bool valid_backend = false;
    #ifdef Heffte_ENABLE_FFTW
    valid_backend = valid_backend or perform_benchmark<backend::fftw>(precision_string, backend_string, "fftw", size_fft, arguments(argc, argv));
    #endif
    valid_backend = valid_backend or perform_benchmark<backend::stock>(precision_string, backend_string, "stock", size_fft, arguments(argc, argv));
    #ifdef Heffte_ENABLE_MKL
    valid_backend = valid_backend or perform_benchmark<backend::mkl>(precision_string, backend_string, "mkl", size_fft, arguments(argc, argv));
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
