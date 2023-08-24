/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFFTE_COS_EXECUTOR_H
#define HEFFFTE_COS_EXECUTOR_H

#include "heffte_pack3d.h"

/*!
 * \ingroup fft3d
 * \addtogroup fft3dr2r Sine and Cosine Transforms
 *
 * HeFFTe now supports the discrete Sine and Cosine transforms,
 * which use real input and output data but the same data layout
 * and MPI communication patterns as the Fourier transform.
 * See heffte::rtransform for more details.
 *
 * The transforms are performed using the standard FFT algorithms
 * with a pre- and post-processing steps.
 */

namespace heffte {

/*!
 * \ingroup fft3dr2r
 * \brief Create a box with larger dimension that will exploit the symmetry for the Sine and Cosine Transforms.
 */
template<typename index>
box3d<index> make_cos_box(box3d<index> const &box){
    std::array<index, 3> high{box.size[0]-1, box.size[1]-1, box.size[2]-1};
    high[box.order[0]] = 4 * box.osize(0) - 1;
    return box3d<index>(std::array<index, 3>{0, 0, 0}, high, box.order);
}

/*!
 * \ingroup fft3dr2r
 * \brief Pre/Post processing for the Cosine transform using the CPU.
 */
struct cpu_cos_pre_pos_processor{
    //! \brief Pre-process in the forward transform.
    template<typename precision>
    static void pre_forward(void*, int length, precision const input[], precision fft_signal[]){
        for(int i = 0; i < length; i++){
            fft_signal[2*i] = 0;
            fft_signal[2*i+1] = input[i];
        }
        fft_signal[2*length] = 0;
        for(int i = 0; i < 2*length; i++){
            fft_signal[4*length-i] = fft_signal[i];
        }
    }
    //! \brief Post-process in the forward transform.
    template<typename precision>
    static void post_forward(void*, int length, std::complex<precision> const fft_result[], precision result[]){
        for(int i = 0; i < length; i++){
            result[i] = std::real(fft_result[i]);
        }
    }
    //! \brief Pre-process in the inverse transform.
    template<typename precision>
    static void pre_backward(void*, int length, precision const input[], std::complex<precision> fft_signal[]){
        for(int i = 0; i < length; i++){
            fft_signal[i] = std::complex<precision>(input[i]);
        }
        fft_signal[length] = 0.0;

        int index = length-1;
        for(int i = length+1; i < 2*length+1; i++){
            fft_signal[i] = std::complex<precision>(-1.0 * input[index]);
            index --;
        }
    }
    //! \brief Post-process in the inverse transform.
    template<typename precision>
    static void post_backward(void*, int length, precision const fft_result[], precision result[]){
        for(int i=0; i<length; i++)
            result[i] = fft_result[2*i + 1];
    }
};

/*!
 * \ingroup fft3dr2r
 * \brief Pre/Post processing for the Sine transform using the CPU.
 */
struct cpu_sin_pre_pos_processor{
    //! \brief Pre-process in the forward transform.
    template<typename precision>
    static void pre_forward(void*, int length, precision const input[], precision fft_signal[]){
        for(int i=0; i<length; i++){
            fft_signal[2*i]   = 0.0;
            fft_signal[2*i+1] = input[i];
        }
        fft_signal[2*length] = 0.;
        for(int i=0; i<length; i++){
            fft_signal[4*length-2*i]  = 0.0;
            fft_signal[4*length-2*i-1]= -input[i];
        }
    }
    //! \brief Post-process in the forward transform.
    template<typename precision>
    static void post_forward(void*, int length, std::complex<precision> const fft_result[], precision result[]){
        for(int i=0; i < length; i++)
            result[i] = -std::imag(fft_result[i+1]);
    }
    //! \brief Pre-process in the inverse transform.
    template<typename precision>
    static void pre_backward(void*, int length, precision const input[], std::complex<precision> fft_signal[]){
        fft_signal[0] = std::complex<precision>(0.0);
        for(int i=0; i < length; i++){
            fft_signal[i+1] = std::complex<precision>(0.0, -input[i]);
        }
        fft_signal[2*length] = std::complex<precision>(0.0);
        for(int i=0; i < length-1; i++){
            fft_signal[length + i + 1] = std::complex<precision>(0.0, -input[length - i - 2]);
        }
    }
    //! \brief Post-process in the inverse transform.
    template<typename precision>
    static void post_backward(void*, int length, precision const fft_result[], precision result[]){
        cpu_cos_pre_pos_processor::post_backward(nullptr, length, fft_result, result);
    }
};

struct cpu_cos1_pre_pos_processor{};
struct cpu_sin1_pre_pos_processor{};

/*!
 * \ingroup fft3dr2r
 * \brief Template algorithm for the Sine and Cosine transforms.
 *
 * \tparam fft_backend_tag indicate the FFT backend to use, e.g., fftw or cufft.
 * \tparam prepost_processor a collection of methods for pre-post processing the data before/after applying the FFT
 */
template<typename fft_backend_tag, typename prepost_processor>
struct real2real_executor : public executor_base{
    //! \brief Construct a plan for batch 1D transforms.
    template<typename index>
    real2real_executor(typename backend::device_instance<typename backend::buffer_traits<fft_backend_tag>::location>::stream_type cstream, box3d<index> const box, int dimension) :
        stream(cstream),
        length(box.osize(0)),
        num_batch(box.osize(1) * box.osize(2)),
        total_size(box.count()),
        fft(make_executor_r2c<fft_backend_tag>(stream, make_cos_box(box), dimension))
    {
        assert(dimension == box.order[0]); // supporting only ordered operations (for now)
    }
    //! \brief Construct a plan for batch 2D transforms, not implemented currently.
    template<typename index>
    real2real_executor(typename backend::device_instance<typename backend::buffer_traits<fft_backend_tag>::location>::stream_type cstream, box3d<index> const, int, int) : stream(cstream)
    { throw std::runtime_error("2D real-to-real transform is not yet implemented!"); }
    //! \brief Construct a plan for a single 3D transform, not implemented currently.
    template<typename index>
    real2real_executor(typename backend::device_instance<typename backend::buffer_traits<fft_backend_tag>::location>::stream_type cstream, box3d<index> const) : stream(cstream)
    { throw std::runtime_error("3D real-to-real transform is not yet implemented!"); }

    //! \brief Forward transform.
    template<typename scalar_type>
    void forward(scalar_type data[], scalar_type workspace[]) const{
        scalar_type* temp = workspace;
        std::complex<scalar_type>* ctemp = align_pntr(reinterpret_cast<std::complex<scalar_type>*>(workspace + fft->box_size() + 1));
        std::complex<scalar_type>* fft_work = (fft->workspace_size() == 0) ? nullptr : ctemp + fft->complex_size();
        for(int i=0; i<num_batch; i++){
            prepost_processor::pre_forward(stream, length, data + i * length, temp + i * 4 * length);
        }
        fft->forward(temp, ctemp, fft_work);
        for(int i=0; i<num_batch; i++)
            prepost_processor::post_forward(stream, length, ctemp + i * (2 * length + 1), data + i * length);
    }
    //! \brief Inverse transform.
    template<typename scalar_type>
    void backward(scalar_type data[], scalar_type workspace[]) const{
        scalar_type* temp = workspace;
        std::complex<scalar_type>* ctemp = align_pntr(reinterpret_cast<std::complex<scalar_type>*>(workspace + fft->box_size() + 1));
        std::complex<scalar_type>* fft_work = (fft->workspace_size() == 0) ? nullptr : ctemp + fft->complex_size();
        for(int i=0; i<num_batch; i++)
            prepost_processor::pre_backward(stream, length, data + i * length, ctemp + i * (2 * length + 1));
        fft->backward(ctemp, temp, fft_work);
        for(int i=0; i<num_batch; i++)
            prepost_processor::post_backward(stream, length, temp + 4 * i * length, data + i * length);
    }

    //! \brief Placeholder for template type consistency, should never be called.
    template<typename precision>
    void forward(precision const[], std::complex<precision>[]) const{
        throw std::runtime_error("Calling cos-transform with real-to-complex data! This should not happen!");
    }
    //! \brief Placeholder for template type consistency, should never be called.
    template<typename precision>
    void backward(std::complex<precision> indata[], precision outdata[]) const{ forward(outdata, indata); }

    //! \brief Returns the size of the box.
    int box_size() const override{ return total_size; }
    //! \brief Returns the size of the box.
    size_t workspace_size() const override{
        return fft->box_size() + 1 + 2 * fft->complex_size() + 2 * fft->workspace_size()
               + ((std::is_same<fft_backend_tag, backend::cufft>::value) ? 1 : 0);
    }
    //! \brief Moves the pointer forward to be aligned to the size of std::complex<scalar_type>, used for CUDA only.
    template<typename scalar_type>
    std::complex<scalar_type>* align_pntr(std::complex<scalar_type> *p) const{
        if (std::is_same<fft_backend_tag, backend::cufft>::value){
            return (reinterpret_cast<size_t>(p) % sizeof(std::complex<scalar_type>) == 0) ? p :
                reinterpret_cast<std::complex<scalar_type>*>(reinterpret_cast<scalar_type*>(p) + 1);
        }else{
            return p;
        }
    }
    //! \brief Forward r2r, single precision.
    virtual void forward(float data[], float *workspace) const override{ forward<float>(data, workspace); }
    //! \brief Forward r2r, double precision.
    virtual void forward(double data[], double *workspace) const override{ forward<double>(data, workspace); }
    //! \brief Backward r2r, single precision.
    virtual void backward(float data[], float *workspace) const override{ backward<float>(data, workspace); }
    //! \brief Backward r2r, double precision.
    virtual void backward(double data[], double *workspace) const override{ backward<double>(data, workspace); }

private:
    typename backend::device_instance<typename backend::buffer_traits<fft_backend_tag>::location>::stream_type stream;

    int length, num_batch, total_size;

    std::unique_ptr<typename one_dim_backend<fft_backend_tag>::executor_r2c> fft;
};


}

#endif
