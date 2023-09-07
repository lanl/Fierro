/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_BACKEND_MKL_H
#define HEFFTE_BACKEND_MKL_H

#include "heffte_r2r_executor.h"

#ifdef Heffte_ENABLE_MKL

#include "mkl_dfti.h"

/*!
 * \ingroup fft3d
 * \addtogroup hefftemkl Backend mkl
 *
 * Wrappers and template specializations related to the MKL backend.
 * Requires CMake option:
 * \code
 *  -D Heffte_ENABLE_MKL=ON
 * \endcode
 */

namespace heffte{

namespace backend{
    /*!
     * \ingroup hefftemkl
     * \brief Indicate that the MKL backend has been enabled.
     */
    template<> struct is_enabled<mkl> : std::true_type{};
    /*!
     * \ingroup hefftemkl
     * \brief Indicate that the MKL Cosine Transform backend has been enabled.
     */
    template<> struct is_enabled<mkl_cos> : std::true_type{};
    /*!
     * \ingroup hefftemkl
     * \brief Indicate that the MKL Sine Transform backend has been enabled.
     */
    template<> struct is_enabled<mkl_sin> : std::true_type{};
    /*!
     * \ingroup hefftemkl
     * \brief Make MKL the default CPU backend.
     */
    template<> struct default_backend<tag::cpu> {
        using type = mkl;
    };
}

/*!
 * \ingroup hefftemkl
 * \brief Recognize the MKL single precision complex type.
 */
template<> struct is_ccomplex<float _Complex> : std::true_type{};
/*!
 * \ingroup hefftemkl
 * \brief Recognize the MKL double precision complex type.
 */
template<> struct is_zcomplex<double _Complex> : std::true_type{};

/*!
 * \ingroup hefftemkl
 * \brief Checks the status of a call to the MKL backend.
 */
inline void check_error(MKL_LONG status, std::string const &function_name){
    if (status != 0){
        throw std::runtime_error(function_name + " failed with status: " + std::to_string(status));
    }
}

/*!
 * \ingroup hefftemkl
 * \brief Base plan for backend::mkl, works only for float and double complex.
 *
 * MKL library uses a unique plan type for forward and backward fft transforms.
 * This class will wrap around such plan and provide RAII style
 * of memory management and simple constructors that take inputs suitable to HeFFTe.
 */
template<typename scalar_type>
struct plan_mkl{
    /*!
     * \brief Constructor, takes inputs identical to MKL FFT descriptors.
     *
     * \param size is the number of entries in a 1-D transform
     * \param howmanyffts is the number of transforms in the batch
     * \param stride is the distance between entries of the same transform
     * \param dist is the distance between the first entries of consecutive sequences
     */
    plan_mkl(int size, int howmanyffts, int stride, int dist) : plan(nullptr){
        static_assert(std::is_same<scalar_type, std::complex<float>>::value
                      or std::is_same<scalar_type, std::complex<double>>::value,
                      "plan_mkl requires std::complex scalar_type with float or double, see plan_mkl_r2c for real types");

        check_error( DftiCreateDescriptor(&plan, (std::is_same<scalar_type, std::complex<float>>::value) ?
                                                  DFTI_SINGLE : DFTI_DOUBLE,
                                          DFTI_COMPLEX, 1, static_cast<MKL_LONG>(size)), "mkl plan create" );
        check_error( DftiSetValue(plan, DFTI_NUMBER_OF_TRANSFORMS, static_cast<MKL_LONG>(howmanyffts)), "mkl set howmany");
        check_error( DftiSetValue(plan, DFTI_PLACEMENT, DFTI_INPLACE), "mkl set in place");
        MKL_LONG lstride[] = {0, static_cast<MKL_LONG>(stride)};
        check_error( DftiSetValue(plan, DFTI_INPUT_STRIDES, lstride), "mkl set istride");
        check_error( DftiSetValue(plan, DFTI_OUTPUT_STRIDES, lstride), "mkl set ostride");
        check_error( DftiSetValue(plan, DFTI_INPUT_DISTANCE, static_cast<MKL_LONG>(dist)), "mkl set idist");
        check_error( DftiSetValue(plan, DFTI_OUTPUT_DISTANCE, static_cast<MKL_LONG>(dist)), "mkl set odist");
        check_error( DftiCommitDescriptor(plan), "mkl commit");
    }
    /*!
     * \brief Constructor, takes inputs identical to MKL FFT descriptors.
     *
     * \param size1 is the number of entries in a 2-D transform, direction 1
     * \param size2 is the number of entries in a 2-D transform, direction 2
     * \param embed is the size of the leading dimensions of the data array
     * \param howmanyffts is the number of transforms in the batch
     * \param dist is the distance between the first entries of consecutive sequences
     */
    plan_mkl(size_t size1, size_t size2, std::array<MKL_LONG, 2> const &embed, size_t howmanyffts, size_t dist) : plan(nullptr){
        static_assert(std::is_same<scalar_type, std::complex<float>>::value
                      or std::is_same<scalar_type, std::complex<double>>::value,
                      "plan_mkl requires std::complex scalar_type with float or double, see plan_mkl_r2c for real types");

        MKL_LONG size[] = {static_cast<MKL_LONG>(size1), static_cast<MKL_LONG>(size2)};
        check_error( DftiCreateDescriptor(&plan, (std::is_same<scalar_type, std::complex<float>>::value) ?
                                                  DFTI_SINGLE : DFTI_DOUBLE,
                                          DFTI_COMPLEX, 2, size), "mkl plan create 2d" );
        check_error( DftiSetValue(plan, DFTI_NUMBER_OF_TRANSFORMS, static_cast<MKL_LONG>(howmanyffts)), "mkl set howmany");
        check_error( DftiSetValue(plan, DFTI_PLACEMENT, DFTI_INPLACE), "mkl set in place");
        MKL_LONG lstride[] = {0, static_cast<MKL_LONG>(embed[0]), static_cast<MKL_LONG>(embed[1])};
        check_error( DftiSetValue(plan, DFTI_INPUT_STRIDES, lstride), "mkl set istride");
        check_error( DftiSetValue(plan, DFTI_OUTPUT_STRIDES, lstride), "mkl set ostride");
        check_error( DftiSetValue(plan, DFTI_INPUT_DISTANCE, static_cast<MKL_LONG>(dist)), "mkl set idist");
        check_error( DftiSetValue(plan, DFTI_OUTPUT_DISTANCE, static_cast<MKL_LONG>(dist)), "mkl set odist");
        check_error( DftiCommitDescriptor(plan), "mkl commit");

    }
    /*!
     * \brief Constructor, takes inputs identical to MKL FFT descriptors.
     *
     * \param size1 is the number of entries in a 3-D transform, direction 1
     * \param size2 is the number of entries in a 3-D transform, direction 2
     * \param size3 is the number of entries in a 3-D transform, direction 3
     */
    plan_mkl(int size1, int size2, int size3){
        MKL_LONG size[] = {static_cast<MKL_LONG>(size3), static_cast<MKL_LONG>(size2), static_cast<MKL_LONG>(size1)};
        check_error( DftiCreateDescriptor(&plan, (std::is_same<scalar_type, std::complex<float>>::value) ?
                                                  DFTI_SINGLE : DFTI_DOUBLE,
                                          DFTI_COMPLEX, 3, size), "mkl plan create 3d" );
        check_error( DftiSetValue(plan, DFTI_NUMBER_OF_TRANSFORMS, 1), "mkl set howmany");
        check_error( DftiSetValue(plan, DFTI_PLACEMENT, DFTI_INPLACE), "mkl set in place");
        check_error( DftiCommitDescriptor(plan), "mkl commit");
    }

    //! \brief Destructor, deletes the plan.
    ~plan_mkl(){ check_error( DftiFreeDescriptor(&plan), "mkl delete descriptor"); }
    //! \brief Custom conversion to the MKL plan.
    operator DFTI_DESCRIPTOR_HANDLE() const{ return plan; }
    //! \brief The MKL opaque structure (pointer to struct).
    DFTI_DESCRIPTOR_HANDLE plan;
};

/*!
 * \ingroup hefftemkl
 * \brief Wrapper around the MKL API.
 *
 * A single class that manages the plans and executions of mkl
 * so that a single API is provided for all backends.
 * The executor operates on a box and performs 1-D FFTs
 * for the given dimension.
 * The class silently manages the plans and buffers needed
 * for the different types.
 * All input and output arrays must have size equal to the box.
 */
class mkl_executor : public executor_base{
public:
    //! \brief Bring forth method that have not been overloaded.
    using executor_base::forward;
    //! \brief Bring forth method that have not been overloaded.
    using executor_base::backward;
    //! \brief Bring forth method that have not been overloaded.
    using executor_base::complex_size;
    //! \brief Constructor, specifies the box and dimension.
    template<typename index>
    mkl_executor(void*, box3d<index> const box, int dimension) :
        size(box.size[dimension]), size2(0),
        howmanyffts(fft1d_get_howmany(box, dimension)),
        stride(fft1d_get_stride(box, dimension)),
        dist((dimension == box.order[0]) ? size : 1),
        blocks((dimension == box.order[1]) ? box.osize(2) : 1),
        block_stride(box.osize(0) * box.osize(1)),
        total_size(box.count()),
        embed({0, 0})
    {}
    //! \brief Merges two FFTs into one.
    template<typename index>
    mkl_executor(void*, box3d<index> const box, int dir1, int dir2) :
        size(box.size[std::min(dir1, dir2)]), size2(box.size[std::max(dir1, dir2)]),
        blocks(1), block_stride(0), total_size(box.count())
    {
        int odir1 = box.find_order(dir1);
        int odir2 = box.find_order(dir2);

        if (std::min(odir1, odir2) == 0 and std::max(odir1, odir2) == 1){
            stride = 1;
            dist = size * size2;
            embed = {static_cast<MKL_LONG>(stride), static_cast<MKL_LONG>(size)};
            howmanyffts = box.size[2];
        }else if (std::min(odir1, odir2) == 1 and std::max(odir1, odir2) == 2){
            stride = box.size[0];
            dist = 1;
            embed = {static_cast<MKL_LONG>(stride), static_cast<MKL_LONG>(size) * static_cast<MKL_LONG>(stride)};
            howmanyffts = box.size[0];
        }else{ // case of directions (0, 2)
            stride = 1;
            dist = size;
            embed = {static_cast<MKL_LONG>(stride), static_cast<MKL_LONG>(box.size[1]) * static_cast<MKL_LONG>(box.size[0])};
            howmanyffts = box.size[1];
        }
    }
    //! \brief Merges three FFTs into one.
    template<typename index>
    mkl_executor(void*, box3d<index> const box) :
        size(box.size[0]), size2(box.size[1]), howmanyffts(box.size[2]),
        stride(0), dist(0),
        blocks(1), block_stride(0),
        total_size(box.count()),
        embed({0, 0})
    {}

    //! \brief Forward fft, float-complex case.
    void forward(std::complex<float> data[], std::complex<float>*) const override{
        make_plan(cplan);
        for(int i=0; i<blocks; i++){
            float _Complex* block_data = reinterpret_cast<float _Complex*>(data + i * block_stride);
            DftiComputeForward(*cplan, block_data);
        }
    }
    //! \brief Backward fft, float-complex case.
    void backward(std::complex<float> data[], std::complex<float>*) const override{
        make_plan(cplan);
        for(int i=0; i<blocks; i++){
            float _Complex* block_data = reinterpret_cast<float _Complex*>(data + i * block_stride);
            DftiComputeBackward(*cplan, block_data);
        }
    }
    //! \brief Forward fft, double-complex case.
    void forward(std::complex<double> data[], std::complex<double>*) const override{
        make_plan(zplan);
        for(int i=0; i<blocks; i++){
            double _Complex* block_data = reinterpret_cast<double _Complex*>(data + i * block_stride);
            DftiComputeForward(*zplan, block_data);
        }
    }
    //! \brief Backward fft, double-complex case.
    void backward(std::complex<double> data[], std::complex<double>*) const override{
        make_plan(zplan);
        for(int i=0; i<blocks; i++){
            double _Complex* block_data = reinterpret_cast<double _Complex*>(data + i * block_stride);
            DftiComputeBackward(*zplan, block_data);
        }
    }

    //! \brief Converts the deal data to complex and performs float-complex forward transform.
    void forward(float const indata[], std::complex<float> outdata[], std::complex<float> *workspace) const override{
        for(int i=0; i<total_size; i++) outdata[i] = std::complex<float>(indata[i]);
        forward(outdata, workspace);
    }
    //! \brief Performs backward float-complex transform and truncates the complex part of the result.
    void backward(std::complex<float> indata[], float outdata[], std::complex<float> *workspace) const override{
        backward(indata, workspace);
        for(int i=0; i<total_size; i++) outdata[i] = std::real(indata[i]);
    }
    //! \brief Converts the deal data to complex and performs double-complex forward transform.
    void forward(double const indata[], std::complex<double> outdata[], std::complex<double> *workspace) const override{
        for(int i=0; i<total_size; i++) outdata[i] = std::complex<double>(indata[i]);
        forward(outdata, workspace);
    }
    //! \brief Performs backward double-complex transform and truncates the complex part of the result.
    void backward(std::complex<double> indata[], double outdata[], std::complex<double> *workspace) const override{
        backward(indata, workspace);
        for(int i=0; i<total_size; i++) outdata[i] = std::real(indata[i]);
    }

    //! \brief Returns the size of the box.
    int box_size() const override{ return total_size; }
    //! \brief Return the size of the needed workspace.
    size_t workspace_size() const override{ return 0; }

private:
    //! \brief Helper template to create the plan.
    template<typename scalar_type>
    void make_plan(std::unique_ptr<plan_mkl<scalar_type>> &plan) const{
        if (not plan){
            if (dist == 0)
                plan = std::unique_ptr<plan_mkl<scalar_type>>(new plan_mkl<scalar_type>(size, size2, howmanyffts));
            else if (size2 == 0)
                plan = std::unique_ptr<plan_mkl<scalar_type>>(new plan_mkl<scalar_type>(size, howmanyffts, stride, dist));
            else
                plan = std::unique_ptr<plan_mkl<scalar_type>>(new plan_mkl<scalar_type>(size, size2, embed, howmanyffts, dist));
        }
    }

    int size, size2, howmanyffts, stride, dist, blocks, block_stride, total_size;
    std::array<MKL_LONG, 2> embed;
    mutable std::unique_ptr<plan_mkl<std::complex<float>>> cplan;
    mutable std::unique_ptr<plan_mkl<std::complex<double>>> zplan;
};

/*!
 * \ingroup hefftemkl
 * \brief Unlike the C2C plan R2C is non-symmetric and it requires that the direction is built into the plan.
 */
template<typename scalar_type, direction dir>
struct plan_mkl_r2c{
    /*!
     * \brief Constructor taking into account the different sizes for the real and complex parts.
     *
     * \param size is the number of entries in a 1-D transform
     * \param howmanyffts is the number of transforms in the batch
     * \param stride is the distance between entries of the same transform
     * \param rdist is the distance between successive 1-D transforms in the real array
     * \param cdist is the distance between successive 1-D transforms in the complex array
     */
    plan_mkl_r2c(int size, int howmanyffts, int stride, int rdist, int cdist) : plan(nullptr){

        static_assert(std::is_same<scalar_type, float>::value or std::is_same<scalar_type, double>::value,
                      "plan_mkl_r2c requires scalar_type with float or double, see plan_mkl for complex types");

        check_error( DftiCreateDescriptor(&plan, (std::is_same<scalar_type,float>::value) ?
                                                  DFTI_SINGLE : DFTI_DOUBLE,
                                          DFTI_REAL, 1, static_cast<MKL_LONG>(size)), "mkl create r2c");
        check_error( DftiSetValue(plan, DFTI_NUMBER_OF_TRANSFORMS, static_cast<MKL_LONG>(howmanyffts)), "mkl set howmany r2c");
        check_error( DftiSetValue(plan, DFTI_PLACEMENT, DFTI_NOT_INPLACE), "mkl set not in place r2c");
        check_error( DftiSetValue(plan, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX), "mkl conj storage cc");
        MKL_LONG lstride[] = {0, static_cast<MKL_LONG>(stride)};
        check_error( DftiSetValue(plan, DFTI_INPUT_STRIDES, lstride), "mkl set istride r2c");
        check_error( DftiSetValue(plan, DFTI_OUTPUT_STRIDES, lstride), "mkl set ostride r2c");
        if (dir == direction::forward){
            check_error( DftiSetValue(plan, DFTI_INPUT_DISTANCE, static_cast<MKL_LONG>(rdist)), "mkl set rdist r2c");
            check_error( DftiSetValue(plan, DFTI_OUTPUT_DISTANCE, static_cast<MKL_LONG>(cdist)), "mkl set cdist r2c");
        }else{
            check_error( DftiSetValue(plan, DFTI_OUTPUT_DISTANCE, static_cast<MKL_LONG>(rdist)), "mkl set back rdist r2c");
            check_error( DftiSetValue(plan, DFTI_INPUT_DISTANCE, static_cast<MKL_LONG>(cdist)), "mkl set back cdist r2c");
        }
        check_error( DftiCommitDescriptor(plan), "mkl commit r2c");
    }


    //! \brief Identical to the float-complex specialization.
    ~plan_mkl_r2c(){ check_error( DftiFreeDescriptor(&plan), "mkl free r2c"); }
    //! \brief Identical to the float-complex specialization.
    operator DFTI_DESCRIPTOR_HANDLE() const{ return plan; }
    //! \brief Identical to the float-complex specialization.
    DFTI_DESCRIPTOR_HANDLE plan;
};

/*!
 * \ingroup hefftemkl
 * \brief Wrapper to mkl API for real-to-complex transform with shortening of the data.
 *
 * Serves the same purpose of heffte::mkl_executor but only real input is accepted
 * and only the unique (non-conjugate) coefficients are computed.
 * All real arrays must have size of real_size() and all complex arrays must have size complex_size().
 */
class mkl_executor_r2c : public executor_base{
public:
    //! \brief Bring forth method that have not been overloaded.
    using executor_base::forward;
    //! \brief Bring forth method that have not been overloaded.
    using executor_base::backward;
    /*!
     * \brief Constructor defines the box and the dimension of reduction.
     *
     * Note that the result sits in the box returned by box.r2c(dimension).
     */
    template<typename index>
    mkl_executor_r2c(void*, box3d<index> const box, int dimension) :
        size(box.size[dimension]),
        howmanyffts(fft1d_get_howmany(box, dimension)),
        stride(fft1d_get_stride(box, dimension)),
        blocks((dimension == box.order[1]) ? box.osize(2) : 1),
        rdist((dimension == box.order[0]) ? size : 1),
        cdist((dimension == box.order[0]) ? size/2 + 1 : 1),
        rblock_stride(box.osize(0) * box.osize(1)),
        cblock_stride(box.osize(0) * (box.osize(1)/2 + 1)),
        rsize(box.count()),
        csize(box.r2c(dimension).count())
    {}

    //! \brief Forward transform, single precision.
    void forward(float const indata[], std::complex<float> outdata[], std::complex<float>*) const override{
        make_plan(sforward);
        for(int i=0; i<blocks; i++){
            float *rdata = const_cast<float*>(indata + i * rblock_stride);
            float _Complex* cdata = reinterpret_cast<float _Complex*>(outdata + i * cblock_stride);
            DftiComputeForward(*sforward, rdata, cdata);
        }
    }
    //! \brief Backward transform, single precision.
    void backward(std::complex<float> indata[], float outdata[], std::complex<float>*) const override{
        make_plan(sbackward);
        for(int i=0; i<blocks; i++){
            float _Complex* cdata = const_cast<float _Complex*>(reinterpret_cast<float _Complex const*>(indata + i * cblock_stride));
            DftiComputeBackward(*sbackward, cdata, outdata + i * rblock_stride);
        }
    }
    //! \brief Forward transform, double precision.
    void forward(double const indata[], std::complex<double> outdata[], std::complex<double>*) const override{
        make_plan(dforward);
        for(int i=0; i<blocks; i++){
            double *rdata = const_cast<double*>(indata + i * rblock_stride);
            double _Complex* cdata = reinterpret_cast<double _Complex*>(outdata + i * cblock_stride);
            DftiComputeForward(*dforward, rdata, cdata);
        }
    }
    //! \brief Backward transform, double precision.
    void backward(std::complex<double> indata[], double outdata[], std::complex<double>*) const override{
        make_plan(dbackward);
        for(int i=0; i<blocks; i++){
            double _Complex* cdata = const_cast<double _Complex*>(reinterpret_cast<double _Complex const*>(indata + i * cblock_stride));
            DftiComputeBackward(*dbackward, cdata, outdata + i * rblock_stride);
        }
    }

    //! \brief Returns the size of the box with real data.
    int box_size() const override{ return rsize; }
    //! \brief Returns the size of the box with complex coefficients.
    int complex_size() const override{ return csize; }
    //! \brief Return the size of the needed workspace.
    size_t workspace_size() const override{ return 0; }

private:
    //! \brief Helper template to initialize the plan.
    template<typename scalar_type, direction dir>
    void make_plan(std::unique_ptr<plan_mkl_r2c<scalar_type, dir>> &plan) const{
        if (!plan) plan = std::unique_ptr<plan_mkl_r2c<scalar_type, dir>>(new plan_mkl_r2c<scalar_type, dir>(size, howmanyffts, stride, rdist, cdist));
    }

    int size, howmanyffts, stride, blocks;
    int rdist, cdist, rblock_stride, cblock_stride, rsize, csize;
    mutable std::unique_ptr<plan_mkl_r2c<float, direction::forward>> sforward;
    mutable std::unique_ptr<plan_mkl_r2c<double, direction::forward>> dforward;
    mutable std::unique_ptr<plan_mkl_r2c<float, direction::backward>> sbackward;
    mutable std::unique_ptr<plan_mkl_r2c<double, direction::backward>> dbackward;
};

/*!
 * \ingroup hefftemkl
 * \brief Helper struct that defines the types and creates instances of one-dimensional executors.
 *
 * The struct is specialized for each backend.
 */
template<> struct one_dim_backend<backend::mkl>{
    //! \brief Defines the complex-to-complex executor.
    using executor = mkl_executor;
    //! \brief Defines the real-to-complex executor.
    using executor_r2c = mkl_executor_r2c;
};
/*!
 * \ingroup hefftemkl
 * \brief Helper struct that defines the types and creates instances of one-dimensional executors.
 *
 * The struct is specialized for each backend.
 */
template<> struct one_dim_backend<backend::mkl_cos>{
    //! \brief Defines the complex-to-complex executor.
    using executor = real2real_executor<backend::mkl, cpu_cos_pre_pos_processor>;
    //! \brief Defines the real-to-complex executor.
    using executor_r2c = void;
};
/*!
 * \ingroup hefftemkl
 * \brief Helper struct that defines the types and creates instances of one-dimensional executors.
 *
 * The struct is specialized for each backend.
 */
template<> struct one_dim_backend<backend::mkl_sin>{
    //! \brief Defines the complex-to-complex executor.
    using executor = real2real_executor<backend::mkl, cpu_sin_pre_pos_processor>;
    //! \brief Defines the real-to-complex executor.
    using executor_r2c = void;
};

/*!
 * \ingroup hefftemkl
 * \brief Sets the default options for the mkl backend.
 */
template<> struct default_plan_options<backend::mkl>{
    //! \brief The reshape operations will not transpose the data.
    static const bool use_reorder = false;
};
/*!
 * \ingroup hefftemkl
 * \brief Sets the default options for the mkl backend.
 */
template<> struct default_plan_options<backend::mkl_cos>{
    //! \brief The reshape operations will not transpose the data.
    static const bool use_reorder = true;
};
/*!
 * \ingroup hefftemkl
 * \brief Sets the default options for the mkl backend.
 */
template<> struct default_plan_options<backend::mkl_sin>{
    //! \brief The reshape operations will not transpose the data.
    static const bool use_reorder = true;
};

}

#endif

#endif   /* HEFFTE_BACKEND_MKL_H */
