/**
 * @class
 * CPU functions of HEFFT
 */
 /*
     -- heFFTe --
        Univ. of Tennessee, Knoxville
        @date
 */

#include "heffte.h"

#include "heffte_c_defines.h"

namespace heffte{

    bool backend_valid_c(int backend){
        if (backend == Heffte_BACKEND_STOCK){
            return true;
        }else
        #ifdef Heffte_ENABLE_FFTW
        if (backend == Heffte_BACKEND_FFTW){
            return true;
        }else
        #endif
        #ifdef Heffte_ENABLE_MKL
        if (backend == Heffte_BACKEND_MKL){
            return true;
        }else
        #endif
        #ifdef Heffte_ENABLE_CUDA
        if (backend == Heffte_BACKEND_CUFFT){
            return true;
        }else
        #endif
        #ifdef Heffte_ENABLE_ROCM
        if (backend == Heffte_BACKEND_ROCFFT){
            return true;
        }else
        #endif
        {
            return false;
        }
    }


    #if defined(Heffte_ENABLE_FFTW)
    using dummy_backend = heffte::backend::fftw;
    #elif defined(Heffte_ENABLE_MKL)
    using dummy_backend = heffte::backend::mkl;
    #elif defined(Heffte_ENABLE_CUDA)
    using dummy_backend = heffte::backend::cufft;
    #elif defined(Heffte_ENABLE_ROCM)
    using dummy_backend = heffte::backend::rocfft;
    #else
    using dummy_backend = heffte::backend::stock;
    #endif

    plan_options options_from_c_options(heffte_plan_options const *options){
        plan_options opts = heffte::default_options<heffte::dummy_backend>();
        opts.use_reorder  = (options->use_reorder != 0);
        opts.use_pencils  = (options->use_pencils != 0);
        switch (options->algorithm){
            case Heffte_RESHAPE_ALGORITHM_ALLTOALLV  : opts.algorithm = reshape_algorithm::alltoallv; break;
            case Heffte_RESHAPE_ALGORITHM_P2P_PLINED : opts.algorithm = reshape_algorithm::p2p_plined; break;
            case Heffte_RESHAPE_ALGORITHM_P2P        : opts.algorithm = reshape_algorithm::p2p; break;
            default:
                throw std::runtime_error("Invalid reshape algorithm in heffte_plan_options");
        };
        opts.use_gpu_aware = (options->use_gpu_aware != 0);
        return opts;
    }

    box3d<> make_box3d_c(int const low[3], int const high[3], int const *order){
        if (order != nullptr){
            return box3d<>({low[0], low[1], low[2]}, {high[0], high[1], high[2]}, {order[0], order[1], order[2]});
        }else{
            return box3d<>({low[0], low[1], low[2]}, {high[0], high[1], high[2]});
        }
    }

    scale get_scaling_c(int s){
        switch(s){
            case Heffte_SCALE_FULL: return scale::full;
            case Heffte_SCALE_SYMMETRIC: return scale::symmetric;
            default:
                return scale::none;
        }
    }

    struct size_inbox_tag{};
    struct size_outbox_tag{};
    struct size_workspace_tag{};
    struct forward_tag{};
    struct backward_tag{};

    template<typename fftclass>
    int call_returnable_c_template(fftclass const *fft, size_inbox_tag const){
        return fft->size_inbox();
    }

    template<typename fftclass>
    int call_returnable_c_template(fftclass const *fft, size_outbox_tag const){
        return fft->size_outbox();
    }

    template<typename fftclass>
    int call_returnable_c_template(fftclass const *fft, size_workspace_tag const){
        return static_cast<int>(fft->size_workspace());
    }

    template<typename fftclass, typename... vars>
    int call_returnable_c_template(fftclass const *fft, forward_tag const, vars... args){
        fft->forward(args...);
        return 0;
    }

    template<typename fftclass, typename... vars>
    int call_returnable_c_template(fftclass const *fft, backward_tag const, vars... args){
        fft->backward(args...);
        return 0;
    }

    template<typename fname, typename... vars>
    int call_returnable_c(heffte_plan const plan, vars... args){
        if (plan->using_r2c){
            if (plan->backend_type == Heffte_BACKEND_STOCK)
                return call_returnable_c_template(reinterpret_cast<heffte::fft3d_r2c<heffte::backend::stock> const*>(plan->fft), fname(), args...);
            #ifdef Heffte_ENABLE_FFTW
            if (plan->backend_type == Heffte_BACKEND_FFTW)
                return call_returnable_c_template(reinterpret_cast<heffte::fft3d_r2c<heffte::backend::fftw> const*>(plan->fft), fname(), args...);
            #endif
            #ifdef Heffte_ENABLE_MKL
            if (plan->backend_type == Heffte_BACKEND_MKL)
                return call_returnable_c_template(reinterpret_cast<heffte::fft3d_r2c<heffte::backend::mkl> const*>(plan->fft), fname(), args...);
            #endif
            #ifdef Heffte_ENABLE_CUDA
            if (plan->backend_type == Heffte_BACKEND_CUFFT)
                return call_returnable_c_template(reinterpret_cast<heffte::fft3d_r2c<heffte::backend::cufft> const*>(plan->fft), fname(), args...);
            #endif
            #ifdef Heffte_ENABLE_ROCM
            if (plan->backend_type == Heffte_BACKEND_ROCFFT)
                return call_returnable_c_template(reinterpret_cast<heffte::fft3d_r2c<heffte::backend::rocfft> const*>(plan->fft), fname(), args...);
            #endif
        }else{
            if (plan->backend_type == Heffte_BACKEND_STOCK)
                return call_returnable_c_template(reinterpret_cast<heffte::fft3d<heffte::backend::stock> const*>(plan->fft), fname(), args...);
            #ifdef Heffte_ENABLE_FFTW
            if (plan->backend_type == Heffte_BACKEND_FFTW)
                return call_returnable_c_template(reinterpret_cast<heffte::fft3d<heffte::backend::fftw> const*>(plan->fft), fname(), args...);
            #endif
            #ifdef Heffte_ENABLE_MKL
            if (plan->backend_type == Heffte_BACKEND_MKL)
                return call_returnable_c_template(reinterpret_cast<heffte::fft3d<heffte::backend::mkl> const*>(plan->fft), fname(), args...);
            #endif
            #ifdef Heffte_ENABLE_CUDA
            if (plan->backend_type == Heffte_BACKEND_CUFFT)
                return call_returnable_c_template(reinterpret_cast<heffte::fft3d<heffte::backend::cufft> const*>(plan->fft), fname(), args...);
            #endif
            #ifdef Heffte_ENABLE_ROCM
            if (plan->backend_type == Heffte_BACKEND_ROCFFT)
                return call_returnable_c_template(reinterpret_cast<heffte::fft3d<heffte::backend::rocfft> const*>(plan->fft), fname(), args...);
            #endif
        }
        return call_returnable_c_template(reinterpret_cast<heffte::fft3d<dummy_backend> const*>(plan->fft), fname(), args...);
    }

    template<typename fname, typename... vars>
    int call_returnable_c_nor2c(heffte_plan const plan, vars... args){
        if (plan->backend_type == Heffte_BACKEND_STOCK)
            return call_returnable_c_template(reinterpret_cast<heffte::fft3d<heffte::backend::stock> const*>(plan->fft), fname(), args...);
        #ifdef Heffte_ENABLE_FFTW
        if (plan->backend_type == Heffte_BACKEND_FFTW)
            return call_returnable_c_template(reinterpret_cast<heffte::fft3d<heffte::backend::fftw> const*>(plan->fft), fname(), args...);
        #endif
        #ifdef Heffte_ENABLE_MKL
        if (plan->backend_type == Heffte_BACKEND_MKL)
            return call_returnable_c_template(reinterpret_cast<heffte::fft3d<heffte::backend::mkl> const*>(plan->fft), fname(), args...);
        #endif
        #ifdef Heffte_ENABLE_CUDA
        if (plan->backend_type == Heffte_BACKEND_CUFFT)
            return call_returnable_c_template(reinterpret_cast<heffte::fft3d<heffte::backend::cufft> const*>(plan->fft), fname(), args...);
        #endif
        #ifdef Heffte_ENABLE_ROCM
        if (plan->backend_type == Heffte_BACKEND_ROCFFT)
            return call_returnable_c_template(reinterpret_cast<heffte::fft3d<heffte::backend::rocfft> const*>(plan->fft), fname(), args...);
        #endif
        return call_returnable_c_template(reinterpret_cast<heffte::fft3d<dummy_backend> const*>(plan->fft), fname(), args...);
    }

}

extern "C"{

int heffte_set_default_options(int backend, heffte_plan_options *options){
    if (not heffte::backend_valid_c(backend)) return 1;
    heffte::plan_options opts = [=]()->heffte::plan_options{
        if (backend == Heffte_BACKEND_STOCK){
            return heffte::default_options<heffte::backend::stock>();
        }
        #ifdef Heffte_ENABLE_FFTW
        if (backend == Heffte_BACKEND_FFTW){
            return heffte::default_options<heffte::backend::fftw>();
        }
        #endif
        #ifdef Heffte_ENABLE_MKL
        if (backend == Heffte_BACKEND_MKL){
            return heffte::default_options<heffte::backend::mkl>();
        }
        #endif
        #ifdef Heffte_ENABLE_CUDA
        if (backend == Heffte_BACKEND_CUFFT){
            return heffte::default_options<heffte::backend::cufft>();
        }
        #endif
        #ifdef Heffte_ENABLE_ROCM
        if (backend == Heffte_BACKEND_ROCFFT){
            return heffte::default_options<heffte::backend::rocfft>();
        }
        #endif
        return heffte::default_options<heffte::dummy_backend>(); // will never happen
    }();
    options->use_reorder = (opts.use_reorder) ? 1 : 0;
    options->algorithm = Heffte_RESHAPE_ALGORITHM_ALLTOALLV;
    options->use_pencils = (opts.use_pencils) ? 1 : 0;
    return Heffte_SUCCESS;
}

int heffte_plan_create(int backend, int const inbox_low[3], int const inbox_high[3], int const *inbox_order,
                       int const outbox_low[3], int const outbox_high[3], int const *outbox_order,
                       MPI_Comm const comm, heffte_plan_options const *options, heffte_plan *plan){
    if (not heffte::backend_valid_c(backend)) return 1;
    *plan = new heffte_fft_plan;
    (*plan)->backend_type = backend;
    (*plan)->using_r2c = 0;

    heffte::box3d<> inbox  = heffte::make_box3d_c(inbox_low, inbox_high, inbox_order);
    heffte::box3d<> outbox = heffte::make_box3d_c(outbox_low, outbox_high, outbox_order);

    heffte_plan_options const *opts = options;
    heffte_plan_options c_opts;
    if (opts == nullptr){
        heffte_set_default_options(backend, &c_opts);
        opts = &c_opts;
    }
    heffte::plan_options cpp_opts = heffte::options_from_c_options(opts);

    try{
    if (backend == Heffte_BACKEND_STOCK){
        (*plan)->fft = reinterpret_cast<void*>(new heffte::fft3d<heffte::backend::stock>(inbox, outbox, comm, cpp_opts));
    }
    #ifdef Heffte_ENABLE_FFTW
    if (backend == Heffte_BACKEND_FFTW){
        (*plan)->fft = reinterpret_cast<void*>(new heffte::fft3d<heffte::backend::fftw>(inbox, outbox, comm, cpp_opts));
    }
    #endif
    #ifdef Heffte_ENABLE_MKL
    if (backend == Heffte_BACKEND_MKL){
        (*plan)->fft = reinterpret_cast<void*>(new heffte::fft3d<heffte::backend::mkl>(inbox, outbox, comm, cpp_opts));
    }
    #endif
    #ifdef Heffte_ENABLE_CUDA
    if (backend == Heffte_BACKEND_CUFFT){
        (*plan)->fft = reinterpret_cast<void*>(new heffte::fft3d<heffte::backend::cufft>(inbox, outbox, comm, cpp_opts));
    }
    #endif
    #ifdef Heffte_ENABLE_ROCM
    if (backend == Heffte_BACKEND_ROCFFT){
        (*plan)->fft = reinterpret_cast<void*>(new heffte::fft3d<heffte::backend::rocfft>(inbox, outbox, comm, cpp_opts));
    }
    #endif

    }catch(std::runtime_error &){
        delete *plan;
        plan = nullptr;
        return 2;
    }

    return Heffte_SUCCESS;
}

int heffte_plan_create_r2c(int backend, int const inbox_low[3], int const inbox_high[3], int const *inbox_order,
                           int const outbox_low[3], int const outbox_high[3], int const *outbox_order,
                           int r2c_direction, MPI_Comm const comm, heffte_plan_options const *options, heffte_plan *plan){

    if (not heffte::backend_valid_c(backend)) return 1;
    *plan = new heffte_fft_plan;
    (*plan)->backend_type = backend;
    (*plan)->using_r2c = 1;

    heffte::box3d<> inbox  = heffte::make_box3d_c(inbox_low, inbox_high, inbox_order);
    heffte::box3d<> outbox = heffte::make_box3d_c(outbox_low, outbox_high, outbox_order);

    heffte_plan_options const *opts = options;
    heffte_plan_options c_opts;
    if (opts == nullptr){
        heffte_set_default_options(backend, &c_opts);
        opts = &c_opts;
    }
    heffte::plan_options cpp_opts = heffte::options_from_c_options(opts);

    try{
    if (backend == Heffte_BACKEND_STOCK){
        (*plan)->fft = reinterpret_cast<void*>(new heffte::fft3d_r2c<heffte::backend::stock>(inbox, outbox, r2c_direction, comm, cpp_opts));
    }
    #ifdef Heffte_ENABLE_FFTW
    if (backend == Heffte_BACKEND_FFTW){
        (*plan)->fft = reinterpret_cast<void*>(new heffte::fft3d_r2c<heffte::backend::fftw>(inbox, outbox, r2c_direction, comm, cpp_opts));
    }
    #endif
    #ifdef Heffte_ENABLE_MKL
    if (backend == Heffte_BACKEND_MKL){
        (*plan)->fft = reinterpret_cast<void*>(new heffte::fft3d_r2c<heffte::backend::mkl>(inbox, outbox, r2c_direction, comm, cpp_opts));
    }
    #endif
    #ifdef Heffte_ENABLE_CUDA
    if (backend == Heffte_BACKEND_CUFFT){
        (*plan)->fft = reinterpret_cast<void*>(new heffte::fft3d_r2c<heffte::backend::cufft>(inbox, outbox, r2c_direction, comm, cpp_opts));
    }
    #endif
    #ifdef Heffte_ENABLE_ROCM
    if (backend == Heffte_BACKEND_ROCFFT){
        (*plan)->fft = reinterpret_cast<void*>(new heffte::fft3d_r2c<heffte::backend::rocfft>(inbox, outbox, r2c_direction, comm, cpp_opts));
    }
    #endif

    }catch(std::runtime_error &){
        delete *plan;
        plan = nullptr;
        return 2;
    }

    return Heffte_SUCCESS;
}

int heffte_plan_destroy(heffte_plan plan){
    if (not heffte::backend_valid_c(plan->backend_type)) return 3;
    if (plan->using_r2c){
        if (plan->backend_type == Heffte_BACKEND_STOCK)
            delete reinterpret_cast<heffte::fft3d_r2c<heffte::backend::stock>*>(plan->fft);
        #ifdef Heffte_ENABLE_FFTW
        if (plan->backend_type == Heffte_BACKEND_FFTW)
            delete reinterpret_cast<heffte::fft3d_r2c<heffte::backend::fftw>*>(plan->fft);
        #endif
        #ifdef Heffte_ENABLE_MKL
        if (plan->backend_type == Heffte_BACKEND_MKL)
            delete reinterpret_cast<heffte::fft3d_r2c<heffte::backend::mkl>*>(plan->fft);
        #endif
        #ifdef Heffte_ENABLE_CUDA
        if (plan->backend_type == Heffte_BACKEND_CUFFT)
            delete reinterpret_cast<heffte::fft3d_r2c<heffte::backend::cufft>*>(plan->fft);
        #endif
        #ifdef Heffte_ENABLE_ROCM
        if (plan->backend_type == Heffte_BACKEND_ROCFFT)
            delete reinterpret_cast<heffte::fft3d_r2c<heffte::backend::rocfft>*>(plan->fft);
        #endif
    }else{
        if (plan->backend_type == Heffte_BACKEND_STOCK)
            delete reinterpret_cast<heffte::fft3d<heffte::backend::stock>*>(plan->fft);
        #ifdef Heffte_ENABLE_FFTW
        if (plan->backend_type == Heffte_BACKEND_FFTW)
            delete reinterpret_cast<heffte::fft3d<heffte::backend::fftw>*>(plan->fft);
        #endif
        #ifdef Heffte_ENABLE_MKL
        if (plan->backend_type == Heffte_BACKEND_MKL)
            delete reinterpret_cast<heffte::fft3d<heffte::backend::mkl>*>(plan->fft);
        #endif
        #ifdef Heffte_ENABLE_CUDA
        if (plan->backend_type == Heffte_BACKEND_CUFFT)
            delete reinterpret_cast<heffte::fft3d<heffte::backend::cufft>*>(plan->fft);
        #endif
        #ifdef Heffte_ENABLE_ROCM
        if (plan->backend_type == Heffte_BACKEND_ROCFFT)
            delete reinterpret_cast<heffte::fft3d<heffte::backend::rocfft>*>(plan->fft);
        #endif
    }
    delete plan;

    return Heffte_SUCCESS;
}

int heffte_size_inbox(heffte_plan const plan){
    return heffte::call_returnable_c<heffte::size_inbox_tag>(plan);
}
int heffte_size_outbox(heffte_plan const plan){
    return heffte::call_returnable_c<heffte::size_outbox_tag>(plan);
}
int heffte_size_workspace(heffte_plan const plan){
    return heffte::call_returnable_c<heffte::size_workspace_tag>(plan);
}
int heffte_get_backend(heffte_plan const plan){
    return plan->backend_type;
}
int heffte_is_r2c(heffte_plan const plan){
    return (plan->using_r2c) ? 1 : 0;
}

void heffte_forward_s2c(heffte_plan const plan, float const *input, void *o, int s){
    std::complex<float> *output = reinterpret_cast<std::complex<float>*>(o);
    heffte::scale scaling = heffte::get_scaling_c(s);
    heffte::call_returnable_c<heffte::forward_tag>(plan, input, output, scaling);
}
void heffte_forward_c2c(heffte_plan const plan, void const *in, void *o, int s){
    std::complex<float> const *input = reinterpret_cast<std::complex<float> const*>(in);
    std::complex<float> *output = reinterpret_cast<std::complex<float>*>(o);
    heffte::scale scaling = heffte::get_scaling_c(s);
    heffte::call_returnable_c_nor2c<heffte::forward_tag>(plan, input, output, scaling);
}
void heffte_forward_d2z(heffte_plan const plan, double const *input, void *o, int s){
    std::complex<double> *output = reinterpret_cast<std::complex<double>*>(o);
    heffte::scale scaling = heffte::get_scaling_c(s);
    heffte::call_returnable_c<heffte::forward_tag>(plan, input, output, scaling);
}
void heffte_forward_z2z(heffte_plan const plan, void const *in, void *o, int s){
    std::complex<double> const *input = reinterpret_cast<std::complex<double> const*>(in);
    std::complex<double> *output = reinterpret_cast<std::complex<double>*>(o);
    heffte::scale scaling = heffte::get_scaling_c(s);
    heffte::call_returnable_c_nor2c<heffte::forward_tag>(plan, input, output, scaling);
}

void heffte_forward_s2c_buffered(heffte_plan const plan, float const *input, void *o, void *w, int s){
    std::complex<float> *output = reinterpret_cast<std::complex<float>*>(o);
    std::complex<float> *workspace = reinterpret_cast<std::complex<float>*>(w);
    heffte::scale scaling = heffte::get_scaling_c(s);
    heffte::call_returnable_c<heffte::forward_tag>(plan, input, output, workspace, scaling);
}
void heffte_forward_c2c_buffered(heffte_plan const plan, void const *in, void *o, void *w, int s){
    std::complex<float> const *input = reinterpret_cast<std::complex<float> const*>(in);
    std::complex<float> *output = reinterpret_cast<std::complex<float>*>(o);
    std::complex<float> *workspace = reinterpret_cast<std::complex<float>*>(w);
    heffte::scale scaling = heffte::get_scaling_c(s);
    heffte::call_returnable_c_nor2c<heffte::forward_tag>(plan, input, output, workspace, scaling);
}
void heffte_forward_d2z_buffered(heffte_plan const plan, double const *input, void *o, void *w, int s){
    std::complex<double> *output = reinterpret_cast<std::complex<double>*>(o);
    std::complex<double> *workspace = reinterpret_cast<std::complex<double>*>(w);
    heffte::scale scaling = heffte::get_scaling_c(s);
    heffte::call_returnable_c<heffte::forward_tag>(plan, input, output, workspace, scaling);
}
void heffte_forward_z2z_buffered(heffte_plan const plan, void const *in, void *o, void *w, int s){
    std::complex<double> const *input = reinterpret_cast<std::complex<double> const*>(in);
    std::complex<double> *output = reinterpret_cast<std::complex<double>*>(o);
    std::complex<double> *workspace = reinterpret_cast<std::complex<double>*>(w);
    heffte::scale scaling = heffte::get_scaling_c(s);
    heffte::call_returnable_c_nor2c<heffte::forward_tag>(plan, input, output, workspace, scaling);
}

void heffte_backward_c2s(heffte_plan const plan, void const *in, float *output, int s){
    std::complex<float> const *input = reinterpret_cast<std::complex<float> const*>(in);
    heffte::scale scaling = heffte::get_scaling_c(s);
    heffte::call_returnable_c<heffte::backward_tag>(plan, input, output, scaling);
}
void heffte_backward_c2c(heffte_plan const plan, void const *in, void *o, int s){
    std::complex<float> const *input = reinterpret_cast<std::complex<float> const*>(in);
    std::complex<float> *output = reinterpret_cast<std::complex<float>*>(o);
    heffte::scale scaling = heffte::get_scaling_c(s);
    heffte::call_returnable_c_nor2c<heffte::backward_tag>(plan, input, output, scaling);
}
void heffte_backward_z2d(heffte_plan const plan, void const *in, double *output, int s){
    std::complex<double> const *input = reinterpret_cast<std::complex<double> const*>(in);
    heffte::scale scaling = heffte::get_scaling_c(s);
    heffte::call_returnable_c<heffte::backward_tag>(plan, input, output, scaling);
}
void heffte_backward_z2z(heffte_plan const plan, void const *in, void *o, int s){
    std::complex<double> const *input = reinterpret_cast<std::complex<double> const*>(in);
    std::complex<double> *output = reinterpret_cast<std::complex<double>*>(o);
    heffte::scale scaling = heffte::get_scaling_c(s);
    heffte::call_returnable_c_nor2c<heffte::backward_tag>(plan, input, output, scaling);
}

void heffte_backward_c2s_buffered(heffte_plan const plan, void const *in, float *output, void *w, int s){
    std::complex<float> const *input = reinterpret_cast<std::complex<float> const*>(in);
    std::complex<float> *workspace = reinterpret_cast<std::complex<float>*>(w);
    heffte::scale scaling = heffte::get_scaling_c(s);
    heffte::call_returnable_c<heffte::backward_tag>(plan, input, output, workspace, scaling);
}
void heffte_backward_c2c_buffered(heffte_plan const plan, void const *in, void *o, void *w, int s){
    std::complex<float> const *input = reinterpret_cast<std::complex<float> const*>(in);
    std::complex<float> *output = reinterpret_cast<std::complex<float>*>(o);
    std::complex<float> *workspace = reinterpret_cast<std::complex<float>*>(w);
    heffte::scale scaling = heffte::get_scaling_c(s);
    heffte::call_returnable_c_nor2c<heffte::backward_tag>(plan, input, output, workspace, scaling);
}
void heffte_backward_z2d_buffered(heffte_plan const plan, void const *in, double *output, void *w, int s){
    std::complex<double> const *input = reinterpret_cast<std::complex<double> const*>(in);
    std::complex<double> *workspace = reinterpret_cast<std::complex<double>*>(w);
    heffte::scale scaling = heffte::get_scaling_c(s);
    heffte::call_returnable_c<heffte::backward_tag>(plan, input, output, workspace, scaling);
}
void heffte_backward_z2z_buffered(heffte_plan const plan, void const *in, void *o, void *w, int s){
    std::complex<double> const *input = reinterpret_cast<std::complex<double> const*>(in);
    std::complex<double> *output = reinterpret_cast<std::complex<double>*>(o);
    std::complex<double> *workspace = reinterpret_cast<std::complex<double>*>(w);
    heffte::scale scaling = heffte::get_scaling_c(s);
    heffte::call_returnable_c_nor2c<heffte::backward_tag>(plan, input, output, workspace, scaling);
}

} // extern "C"
