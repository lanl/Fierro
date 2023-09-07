/** @class */
/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/
#ifndef TEST_COMMON_H
#define TEST_COMMON_H

#include "heffte.h"

#ifdef Heffte_ENABLE_CUDA
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <cufft.h>
#endif

#define tassert(_result_)          \
    if (!(_result_)){              \
        heffte_test_pass = false;  \
        heffte_all_tests = false;  \
        throw std::runtime_error( std::string("mpi rank = ") + std::to_string(heffte::mpi::comm_rank(MPI_COMM_WORLD)) + "  test " \
                                  + heffte_test_name + " in file: " + __FILE__ + " line: " + std::to_string(__LINE__) );          \
    }

#define sassert(_result_)          \
    if (!(_result_)){              \
        heffte_test_pass = false;  \
        heffte_all_tests = false;  \
        throw std::runtime_error( std::string("  test ") \
                                  + heffte_test_name + " in file: " + __FILE__ + " line: " + std::to_string(__LINE__) );          \
    }

using namespace heffte;

using std::cout;
using std::cerr;
using std::endl;
using std::setw;

std::string heffte_test_name;   // the name of the currently running test
bool heffte_test_pass  = true;  // helps in reporting whether the last test passed
bool heffte_all_tests  = true;  // reports total result of all tests

constexpr int pad_type  = 10;
constexpr int pad_large = 50;
constexpr int pad_pass  = 18;
constexpr int pad_all = pad_type + pad_large + pad_pass + 2;

struct using_mpi{};
struct using_nompi{};

template<typename mpi_tag = using_mpi>
struct all_tests{
    all_tests(char const *cname) : name(cname), separator(pad_all, '-'){
        if (std::is_same<mpi_tag, using_nompi>::value or heffte::mpi::comm_rank(MPI_COMM_WORLD) == 0){
            int const pad = pad_all / 2 + name.size() / 2;
            cout << "\n" << separator << "\n";
            cout << setw(pad) << name << "\n";
            cout << separator << "\n\n";
        }
    }
    ~all_tests(){
        if (std::is_same<mpi_tag, using_nompi>::value or heffte::mpi::comm_rank(MPI_COMM_WORLD) == 0){
            int const pad = pad_all / 2 + name.size() / 2 + 3;
            cout << "\n" << separator << "\n";
            cout << setw(pad) << name  + "  " + ((heffte_all_tests) ? "pass" : "fail") << "\n";
            cout << separator << "\n\n";
        }
    }
    std::string name;
    std::string separator;
};

template<typename scalar_variant> std::string get_variant(){ return ""; }
template<> std::string get_variant<float>(){ return "float"; }
template<> std::string get_variant<double>(){ return "double"; }
template<> std::string get_variant<std::complex<float>>(){ return "ccomplex"; }
template<> std::string get_variant<std::complex<double>>(){ return "zcomplex"; }

struct using_alltoall{};
struct using_pointtopoint{};
template<reshape_algorithm variant> std::string get_description(){ return ""; }
template<> std::string get_description<reshape_algorithm::alltoallv>(){ return "heffte::reshape3d_alltoallv"; }
template<> std::string get_description<reshape_algorithm::alltoall>(){ return "heffte::reshape3d_alltoall"; }
template<> std::string get_description<reshape_algorithm::p2p>(){ return "heffte::reshape3d_pointtopoint"; }
template<> std::string get_description<reshape_algorithm::p2p_plined>(){ return "heffte::reshape3d_p2p (plined)"; }

template<typename scalar_variant = int, typename mpi_tag = using_mpi, typename backend_tag = void>
struct current_test{
    current_test(std::string const &name, MPI_Comm const comm) : test_comm(comm){
        static_assert(std::is_same<mpi_tag, using_mpi>::value, "current_test cannot take a comm when using nompi mode");
        heffte_test_name = name;
        heffte_test_pass = true;
        if (std::is_same<mpi_tag, using_mpi>::value) MPI_Barrier(test_comm);
    };
    current_test(std::string const &name) : test_comm(MPI_COMM_NULL){
        static_assert(std::is_same<mpi_tag, using_nompi>::value, "current_test requires a comm when working in mpi mode");
        heffte_test_name = name;
        heffte_test_pass = true;
    };
    ~current_test(){
        if (std::is_same<mpi_tag, using_nompi>::value or heffte::mpi::comm_rank(MPI_COMM_WORLD) == 0){
            cout << setw(pad_type)  << get_variant<scalar_variant>();
            if (std::is_same<backend_tag, void>::value){
                cout << setw(pad_large) << heffte_test_name;
            }else{
                 cout << setw(pad_large) << heffte_test_name + "<" + heffte::backend::name<backend_tag>() + ">";
            }
            cout << setw(pad_pass)  << ((heffte_test_pass) ? "pass" : "fail") << endl;
        }
        if (std::is_same<mpi_tag, using_mpi>::value) MPI_Barrier(test_comm);
    };
    MPI_Comm const test_comm;
};

template<typename T>
inline bool match(std::vector<T> const &a, std::vector<T> const &b){
    if (a.size() != b.size()) return false;
    for(size_t i=0; i<a.size(); i++)
        if (a[i] != b[i]) return false;
    return true;
}
template<typename T>
inline bool match_verbose(std::vector<T> const &a, std::vector<T> const &b){
    if (a.size() != b.size()) return false;
    for(size_t i=0; i<a.size(); i++){
        if (a[i] != b[i]){
            cout << " mismatch in entry i = " << i << "  with a[i] = " << a[i] << " and b[i] = " << b[i] << endl;
            return false;
        }
    }
    return true;
}

template<typename T> struct precision{};
template<> struct precision<float>{ static constexpr float tolerance = 5.E-4; };
template<> struct precision<double>{ static constexpr double tolerance = 1.E-11; };
template<> struct precision<std::complex<float>>{ static constexpr float tolerance = 5.E-4; };
template<> struct precision<std::complex<double>>{ static constexpr double tolerance = 1.E-11; };

template<typename T>
inline bool approx(std::vector<T> const &a, std::vector<T> const &b, double correction = 1.0){
    if (a.size() != b.size()) return false;
    for(size_t i=0; i<a.size(); i++)
        if (std::abs(a[i] - b[i]) * correction > precision<T>::tolerance){
            cout << " error magnitude: " << std::abs(a[i] - b[i]) * correction << endl;
            return false;
        }
    return true;
}

template<typename backend_tag, typename = void>
struct test_traits{
    template<typename T> using container = typename backend::buffer_traits<backend_tag>::template container<T>;
    template<typename T>
    static container<T> load(std::vector<T> const &x){ return x; }
    template<typename T>
    static std::vector<T> unload(container<T> const &x){ return x; }
};
template<>
struct test_traits<tag::cpu, void>{
    template<typename T> using container = std::vector<T>;
    template<typename T>
    static container<T> load(std::vector<T> const &x){ return x; }
    template<typename T>
    static std::vector<T> unload(container<T> const &x){ return x; }
};

template<typename backend_tag> void* make_stream(backend_tag){ return nullptr; } // CPU case
void sync_stream(void*){}
void free_stream(void*){}

#ifdef Heffte_ENABLE_FFTW
using cpu_backend = heffte::backend::fftw;
#else
#ifdef Heffte_ENABLE_MKL
using cpu_backend = heffte::backend::mkl;
#else
using cpu_backend = heffte::backend::stock;
#endif
#endif

template<typename location_tag>
struct default_tag{
    using backend_tag = cpu_backend;
};

#ifdef Heffte_ENABLE_CUDA
using gpu_backend = heffte::backend::cufft;

cudaStream_t make_stream(backend::cufft){
    cudaStream_t result;
    cudaStreamCreateWithFlags(&result, cudaStreamNonBlocking);
    return result;
}
void sync_stream(cudaStream_t stream){ cudaStreamSynchronize(stream); }
void free_stream(cudaStream_t stream){ cudaStreamDestroy(stream); }
#endif
#ifdef Heffte_ENABLE_ROCM
using gpu_backend = heffte::backend::rocfft;

hipStream_t make_stream(backend::rocfft){
    hipStream_t result;
    hipStreamCreateWithFlags(&result, hipStreamNonBlocking);
    return result;
}
void sync_stream(hipStream_t stream){ hipStreamSynchronize(stream); }
void free_stream(hipStream_t stream){ hipStreamDestroy(stream); }
#endif
#ifdef Heffte_ENABLE_ONEAPI
using gpu_backend = heffte::backend::onemkl;

sycl::queue make_stream(backend::onemkl){
    return heffte::oapi::make_sycl_queue();
}
void sync_stream(sycl::queue &stream){ stream.wait(); }
void free_stream(sycl::queue&){}
#endif
#ifdef Heffte_ENABLE_GPU
template<> struct default_tag<tag::gpu>{
    using backend_tag = gpu_backend;
};
template<typename T>
inline bool match(heffte::gpu::vector<T> const &a, std::vector<T> const &b){
    return match(heffte::gpu::transfer::unload(a), b);
}
template<typename T>
inline bool approx(std::vector<T> const &a, heffte::gpu::vector<T> const &b, double correction = 1.0){
    return approx(a, heffte::gpu::transfer().unload(b), correction);
}
template<typename T>
inline bool approx(heffte::gpu::vector<T> const &a, std::vector<T> const &b, double correction = 1.0){
    return approx(heffte::gpu::transfer().unload(a), b, correction);
}
template<typename T>
inline bool approx(heffte::gpu::vector<T> const &a, heffte::gpu::vector<T> const &b, double correction = 1.0){
    return approx(a, heffte::gpu::transfer::unload(b), correction);
}
template<typename backend_tag>
struct test_traits<backend_tag, typename std::enable_if<backend::uses_gpu<backend_tag>::value, void>::type>{
    template<typename T> using container = gpu::vector<T>;
    template<typename T>
    static container<T> load(std::vector<T> const &x){ return gpu::transfer().load(x); }
    template<typename T>
    static std::vector<T> unload(container<T> const &x){ return gpu::transfer().unload(x); }
};
template<>
struct test_traits<tag::gpu, void>{
    template<typename T> using container = gpu::vector<T>;
    template<typename T>
    static container<T> load(std::vector<T> const &x){ return gpu::transfer().load(x); }
    template<typename T>
    static std::vector<T> unload(container<T> const &x){ return gpu::transfer().unload(x); }
};
#endif

//! \brief Converts a set of c-style strings into a single collection (deque) of c++ strings.
inline std::deque<std::string> arguments(int argc, char *argv[]){
    std::deque<std::string> args;
    for(int i=0; i<argc; i++){
        args.push_back(std::string(argv[i]));
    }
    return args;
}

//! \brief Returns the integer for the subcomm option.
int get_subcomm(std::deque<std::string> const &args){
    auto iopt = args.begin();
    while(iopt != args.end()){
        if (*iopt == "-subcomm"){
            if (++iopt != args.end()){
                return std::stoi(*iopt);
            }else{
                throw std::runtime_error("-subcomm must be followed by an integer");
            }
        }
        iopt++;
    }
    return -1;
}

//! \brief Sets the default options for the backend and then modifies those according to the passed arguments.
template<typename backend_tag>
heffte::plan_options args_to_options(std::deque<std::string> const &args){
    heffte::plan_options options = heffte::default_options<backend_tag>();
    for(auto const &s : args){
        if (s == "-reorder"){
            options.use_reorder = true;
        }else if (s == "-no-reorder"){
            options.use_reorder = false;
        }else if (s == "-a2a"){
            options.algorithm = reshape_algorithm::alltoall;
        }else if (s == "-a2av"){
            options.algorithm = reshape_algorithm::alltoallv;
        }else if (s == "-p2p"){
            options.algorithm = reshape_algorithm::p2p;
        }else if (s == "-p2p_pl"){
            options.algorithm = reshape_algorithm::p2p_plined;
        }else if (s == "-pencils"){
            options.use_pencils = true;
        }else if (s == "-slabs"){
            options.use_pencils = false;
        }else if (s == "-no-gpu-aware"){
            options.use_gpu_aware = false;
        }
    }
    int subcomm = get_subcomm(args);
    if (subcomm != -1) options.use_num_subranks(subcomm);
    return options;
}

template<typename backend_tag>
std::vector<heffte::plan_options> make_all_options(){
    std::vector<heffte::plan_options> result;
    for(int shape = 0; shape < 2; shape++){
        for(int reorder = 0; reorder < 2; reorder++){
            for(reshape_algorithm alg : std::array<reshape_algorithm, 2>{reshape_algorithm::alltoallv, reshape_algorithm::p2p}){
                heffte::plan_options options = default_options<backend_tag>();
                options.use_pencils = (shape == 0);
                options.use_reorder = (reorder == 0);
                options.algorithm = alg;
                result.push_back(options);
            }
        }
    }
    return result;
}

template<typename backend_tag>
std::vector<heffte::plan_options> make_all_options2(){
    std::vector<heffte::plan_options> result;
    for(int shape = 0; shape < 2; shape++){
        for(int reorder = 0; reorder < 2; reorder++){
            for(reshape_algorithm alg : std::array<reshape_algorithm, 3>{reshape_algorithm::alltoall, reshape_algorithm::alltoallv, reshape_algorithm::p2p}){
                heffte::plan_options options = default_options<backend_tag>();
                options.use_pencils = (shape == 0);
                options.use_reorder = (reorder == 0);
                options.algorithm = alg;
                result.push_back(options);
            }
        }
    }
    return result;
}

//! \brief If input and output grid of processors are pencils, useful for comparison with other libraries
bool has_option(std::deque<std::string> const &args, std::string const &opt){
    for(auto &s : args)
        if (s == opt)
            return true;
    return false;
}
//! \brief Takes the three arguments after \b opt and converts them to an array of ints, throws runtime_error if no arguments or cannot convert.
std::array<int, 3> get_grid(std::deque<std::string> const &args, std::string const &opt){
    auto iopt = args.begin();
    while(iopt != args.end()){
        if (*iopt == opt){ // found the argument, take the next three entries
            std::array<int, 3> result = {0, 0, 0};
            for(size_t i=0; i<3; i++){
                if (++iopt != args.end())
                    result[i] = std::stoi(*iopt);
                else
                    throw std::runtime_error(opt + " must be followed by three integers indicating the grid dimensions");
            }
            return result;
        }
        iopt++;
    }
    throw std::runtime_error(opt + " not found");
}

int get_int_arg(std::string const &name, std::deque<std::string> const &args, int default_value = -1){
    auto iopt = args.begin();
    while(iopt != args.end()){
        if (*iopt == name){
            if (++iopt != args.end()){
                return std::stoi(*iopt);
            }else{
                throw std::runtime_error(name + " must be followed by an integer");
            }
        }
        iopt++;
    }
    return default_value;
}

int nruns(std::deque<std::string> const &args){
    for(auto &s : args)
        if (s == "-n10")
            return 10;
        else if (s == "-n50")
            return 50;
    return 5;
}

#endif
