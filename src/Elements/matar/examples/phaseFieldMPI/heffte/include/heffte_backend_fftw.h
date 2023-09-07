/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_BACKEND_FFTW_H
#define HEFFTE_BACKEND_FFTW_H

#include "heffte_r2r_executor.h"

#ifdef Heffte_ENABLE_FFTW

#include "fftw3.h"

/*!
 * \ingroup fft3d
 * \addtogroup hefftefftw Backend fftw3
 *
 * Wrappers and template specializations related to the FFTW backend.
 * Requires CMake option:
 * \code
 *  -D Heffte_ENABLE_FFTW=ON
 * \endcode
 */

namespace heffte{

namespace backend{
    /*!
     * \ingroup hefftefftw
     * \brief Indicate that the FFTW backend has been enabled.
     */
    template<> struct is_enabled<fftw> : std::true_type{};

    /*!
     * \ingroup hefftefftw
     * \brief Indicate that the cos() transform using the FFTW backend has been enabled.
     */
    template<> struct is_enabled<fftw_cos> : std::true_type{};
    /*!
     * \ingroup hefftefftw
     * \brief Indicate that the cos() transform using the FFTW backend has been enabled.
     */
    template<> struct is_enabled<fftw_sin> : std::true_type{};
    /*!
     * \ingroup hefftefftw
     * \brief Indicate that the cos() transform type II using the FFTW backend has been enabled.
     */
    template<> struct is_enabled<fftw_cos1> : std::true_type{};
    /*!
     * \ingroup hefftefftw
     * \brief Indicate that the cos() transform type II using the FFTW backend has been enabled.
     */
    template<> struct is_enabled<fftw_sin1> : std::true_type{};

// Specialization is not necessary since the default behavior assumes CPU parameters.
//     template<>
//     struct buffer_traits<fftw>{
//         using location = tag::cpu;
//         template<typename T> using container = std::vector<T>;
//     };
#ifndef Heffte_ENABLE_MKL
    /*!
     * \ingroup hefftefftw
     * \brief Make FFTW the default CPU backend, if MKL is not enabled.
     */
    template<> struct default_backend<tag::cpu> {
        using type = fftw;
    };
#endif
}

/*!
 * \ingroup hefftefftw
 * \brief Recognize the FFTW single precision complex type.
 */
template<> struct is_ccomplex<fftwf_complex> : std::true_type{};
/*!
 * \ingroup hefftefftw
 * \brief Recognize the FFTW double precision complex type.
 */
template<> struct is_zcomplex<fftw_complex> : std::true_type{};

/*!
 * \ingroup hefftefftw
 * \brief Base plan for fftw, using only the specialization for float and double complex.
 *
 * FFTW3 library uses plans for forward and backward fft transforms.
 * The specializations to this struct will wrap around such plans and provide RAII style
 * of memory management and simple constructors that take inputs suitable to HeFFTe.
 */
template<typename, direction> struct plan_fftw{};

/*!
 * \ingroup hefftefftw
 * \brief Plan for the single precision complex transform.
 *
 * \tparam dir indicates a forward or backward transform
 */
template<direction dir>
struct plan_fftw<std::complex<float>, dir>{
    /*!
     * \brief Constructor, takes inputs identical to fftwf_plan_many_dft().
     *
     * \param size is the number of entries in a 1-D transform
     * \param howmanyffts is the number of transforms in the batch
     * \param stride is the distance between entries of the same transform
     * \param dist is the distance between the first entries of consecutive sequences
     */
    plan_fftw(int size, int howmanyffts, int stride, int dist) :
        plan(fftwf_plan_many_dft(1, &size, howmanyffts, nullptr, nullptr, stride, dist,
                                                    nullptr, nullptr, stride, dist,
                                                    (dir == direction::forward) ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE
                                ))
        {}
    /*!
     * \brief Constructor, takes inputs identical to fftwf_plan_many_dft().
     *
     * \param size1 is the number of entries in a 2-D transform, dimension 1
     * \param size2 is the number of entries in a 2-D transform, dimension 2
     * \param embed is the size of the leading dimensions of the data array
     * \param howmanyffts is the number of transforms in the batch
     * \param stride is the distance between entries of the same transform
     * \param dist is the distance between the first entries of consecutive sequences
     */
    plan_fftw(int size1, int size2, std::array<int, 2> const &embed, int howmanyffts, int stride, int dist){
        std::array<int, 2> size = {size2, size1};

        if (embed[0] == 0 and embed[1] == 0){
            plan = fftwf_plan_many_dft(2, size.data(), howmanyffts, nullptr, nullptr, stride, dist,
                                                    nullptr, nullptr, stride, dist,
                                                    (dir == direction::forward) ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE
                                      );
        }else{
            plan = fftwf_plan_many_dft(2, size.data(), howmanyffts, nullptr, embed.data(), stride, dist,
                                                    nullptr, embed.data(), stride, dist,
                                                    (dir == direction::forward) ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE
                                      );
        }
    }
    //! \brief Identical to the float-complex specialization.
    plan_fftw(int size1, int size2, int size3){
        std::array<int, 3> size = {size3, size2, size1};
        plan = fftwf_plan_many_dft(3, size.data(), 1, nullptr, nullptr, 1, 1, nullptr, nullptr, 1, 1,
                                   (dir == direction::forward) ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
    }
    //! \brief Destructor, deletes the plan.
    ~plan_fftw(){ fftwf_destroy_plan(plan); }
    //! \brief Custom conversion to the FFTW3 plan.
    operator fftwf_plan() const{ return plan; }
    //! \brief The FFTW3 opaque structure (pointer to struct).
    fftwf_plan plan;
};
/*!
 * \ingroup hefftefftw
 * \brief Specialization for double complex.
 */
template<direction dir>
struct plan_fftw<std::complex<double>, dir>{
    //! \brief Identical to the float-complex specialization.
    plan_fftw(int size, int howmanyffts, int stride, int dist) :
        plan(fftw_plan_many_dft(1, &size, howmanyffts, nullptr, nullptr, stride, dist,
                                                   nullptr, nullptr, stride, dist,
                                                   (dir == direction::forward) ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE
                               ))
        {}
    //! \brief Identical to the float-complex specialization.
    plan_fftw(int size1, int size2, std::array<int, 2> const &embed, int howmanyffts, int stride, int dist){
        std::array<int, 2> size = {size2, size1};

        if (embed[0] == 0 and embed[1] == 0){
            plan = fftw_plan_many_dft(2, size.data(), howmanyffts, nullptr, nullptr, stride, dist,
                                                      nullptr, nullptr, stride, dist,
                                                      (dir == direction::forward) ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE
                                     );
        }else{
            plan = fftw_plan_many_dft(2, size.data(), howmanyffts, nullptr, embed.data(), stride, dist,
                                                      nullptr, embed.data(), stride, dist,
                                                      (dir == direction::forward) ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE
                                     );
        }
    }
    //! \brief Identical to the float-complex specialization.
    plan_fftw(int size1, int size2, int size3){
        std::array<int, 3> size = {size3, size2, size1};
        plan = fftw_plan_many_dft(3, size.data(), 1, nullptr, nullptr, 1, 1, nullptr, nullptr, 1, 1,
                                  (dir == direction::forward) ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
    }
    //! \brief Identical to the float-complex specialization.
    ~plan_fftw(){ fftw_destroy_plan(plan); }
    //! \brief Identical to the float-complex specialization.
    operator fftw_plan() const{ return plan; }
    //! \brief Identical to the float-complex specialization.
    fftw_plan plan;
};

/*!
 * \ingroup hefftefftw
 * \brief Wrapper around the FFTW3 API.
 *
 * A single class that manages the plans and executions of fftw3
 * so that a single API is provided for all backends.
 * The executor operates on a box and performs 1-D FFTs
 * for the given dimension.
 * The class silently manages the plans and buffers needed
 * for the different types.
 * All input and output arrays must have size equal to the box.
 */
class fftw_executor : public executor_base{
public:
    //! \brief Bring forth method that have not been overloaded.
    using executor_base::forward;
    //! \brief Bring forth method that have not been overloaded.
    using executor_base::backward;
    //! \brief Bring forth method that have not been overloaded.
    using executor_base::complex_size;
    //! \brief Constructor, specifies the box and dimension.
    template<typename index>
    fftw_executor(void*, box3d<index> const box, int dimension) :
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
    fftw_executor(void*, box3d<index> const box, int dir1, int dir2) :
        size(box.size[std::min(dir1, dir2)]), size2(box.size[std::max(dir1, dir2)]),
        blocks(1), block_stride(0), total_size(box.count()), embed({0, 0})
    {
        int odir1 = box.find_order(dir1);
        int odir2 = box.find_order(dir2);

        if (std::min(odir1, odir2) == 0 and std::max(odir1, odir2) == 1){
            stride = 1;
            dist = size * size2;
            howmanyffts = box.size[2];
        }else if (std::min(odir1, odir2) == 1 and std::max(odir1, odir2) == 2){
            stride = box.size[0];
            dist = 1;
            howmanyffts = box.size[0];
        }else{ // case of directions (0, 2)
            stride = 1;
            dist = size;
            embed = {static_cast<int>(box.size[2]), static_cast<int>(box.size[1] * box.size[0])};
            howmanyffts = box.size[1];
        }
    }
    //! \brief Merges three FFTs into one.
    template<typename index>
    fftw_executor(void*, box3d<index> const box) :
        size(box.size[0]), size2(box.size[1]), howmanyffts(box.size[2]),
        stride(0), dist(0),
        blocks(1), block_stride(0),
        total_size(box.count()),
        embed({0, 0})
    {}

    //! \brief Forward fft, float-complex case.
    void forward(std::complex<float> data[], std::complex<float>*) const override{
        make_plan(cforward);
        for(int i=0; i<blocks; i++){
            fftwf_complex* block_data = reinterpret_cast<fftwf_complex*>(data + i * block_stride);
            fftwf_execute_dft(*cforward, block_data, block_data);
        }
    }
    //! \brief Backward fft, float-complex case.
    void backward(std::complex<float> data[], std::complex<float>*) const override{
        make_plan(cbackward);
        for(int i=0; i<blocks; i++){
            fftwf_complex* block_data = reinterpret_cast<fftwf_complex*>(data + i * block_stride);
            fftwf_execute_dft(*cbackward, block_data, block_data);
        }
    }
    //! \brief Forward fft, double-complex case.
    void forward(std::complex<double> data[], std::complex<double>*) const override{
        make_plan(zforward);
        for(int i=0; i<blocks; i++){
            fftw_complex* block_data = reinterpret_cast<fftw_complex*>(data + i * block_stride);
            fftw_execute_dft(*zforward, block_data, block_data);
        }
    }
    //! \brief Backward fft, double-complex case.
    void backward(std::complex<double> data[], std::complex<double>*) const override{
        make_plan(zbackward);
        for(int i=0; i<blocks; i++){
            fftw_complex* block_data = reinterpret_cast<fftw_complex*>(data + i * block_stride);
            fftw_execute_dft(*zbackward, block_data, block_data);
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
    template<typename scalar_type, direction dir>
    void make_plan(std::unique_ptr<plan_fftw<scalar_type, dir>> &plan) const{
        if (not plan){
            if (dist == 0)
                plan = std::unique_ptr<plan_fftw<scalar_type, dir>>(new plan_fftw<scalar_type, dir>(size, size2, howmanyffts));
            else if (size2 == 0)
                plan = std::unique_ptr<plan_fftw<scalar_type, dir>>(new plan_fftw<scalar_type, dir>(size, howmanyffts, stride, dist));
            else
                plan = std::unique_ptr<plan_fftw<scalar_type, dir>>(new plan_fftw<scalar_type, dir>(size, size2, embed, howmanyffts, stride, dist));
        }
    }

    int size, size2, howmanyffts, stride, dist, blocks, block_stride, total_size;
    std::array<int, 2> embed;
    mutable std::unique_ptr<plan_fftw<std::complex<float>, direction::forward>> cforward;
    mutable std::unique_ptr<plan_fftw<std::complex<float>, direction::backward>> cbackward;
    mutable std::unique_ptr<plan_fftw<std::complex<double>, direction::forward>> zforward;
    mutable std::unique_ptr<plan_fftw<std::complex<double>, direction::backward>> zbackward;
};

/*!
 * \ingroup hefftefftw
 * \brief Specialization for r2c single precision.
 */
template<direction dir>
struct plan_fftw<float, dir>{
    /*!
     * \brief Constructor taking into account the different sizes for the real and complex parts.
     *
     * \param size is the number of entries in a 1-D transform
     * \param howmanyffts is the number of transforms in the batch
     * \param stride is the distance between entries of the same transform
     * \param rdist is the distance between the first entries of consecutive sequences in the real sequences
     * \param cdist is the distance between the first entries of consecutive sequences in the complex sequences
     */
    plan_fftw(int size, int howmanyffts, int stride, int rdist, int cdist) :
        plan((dir == direction::forward) ?
             fftwf_plan_many_dft_r2c(1, &size, howmanyffts, nullptr, nullptr, stride, rdist,
                                                   nullptr, nullptr, stride, cdist,
                                                   FFTW_ESTIMATE
                                   ) :
             fftwf_plan_many_dft_c2r(1, &size, howmanyffts, nullptr, nullptr, stride, cdist,
                                                   nullptr, nullptr, stride, rdist,
                                                   FFTW_ESTIMATE
                                   ))
        {}
    //! \brief Identical to the float-complex specialization.
    ~plan_fftw(){ fftwf_destroy_plan(plan); }
    //! \brief Identical to the float-complex specialization.
    operator fftwf_plan() const{ return plan; }
    //! \brief Identical to the float-complex specialization.
    fftwf_plan plan;
};
/*!
 * \ingroup hefftefftw
 * \brief Specialization for r2c double precision.
 */
template<direction dir>
struct plan_fftw<double, dir>{
    //! \brief Identical to the float-complex specialization.
    plan_fftw(int size, int howmanyffts, int stride, int rdist, int cdist) :
        plan((dir == direction::forward) ?
             fftw_plan_many_dft_r2c(1, &size, howmanyffts, nullptr, nullptr, stride, rdist,
                                                   nullptr, nullptr, stride, cdist,
                                                   FFTW_ESTIMATE
                                   ) :
             fftw_plan_many_dft_c2r(1, &size, howmanyffts, nullptr, nullptr, stride, cdist,
                                                   nullptr, nullptr, stride, rdist,
                                                   FFTW_ESTIMATE
                                   ))
        {}
    //! \brief Identical to the float-complex specialization.
    ~plan_fftw(){ fftw_destroy_plan(plan); }
    //! \brief Identical to the float-complex specialization.
    operator fftw_plan() const{ return plan; }
    //! \brief Identical to the float-complex specialization.
    fftw_plan plan;
};

/*!
 * \ingroup hefftefftw
 * \brief Wrapper to fftw3 API for real-to-complex transform with shortening of the data.
 *
 * Serves the same purpose of heffte::fftw_executor but only real input is accepted
 * and only the unique (non-conjugate) coefficients are computed.
 * All real arrays must have size of real_size() and all complex arrays must have size complex_size().
 */
class fftw_executor_r2c : public executor_base{
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
    fftw_executor_r2c(void*, box3d<index> const box, int dimension) :
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
            fftwf_complex* cdata = reinterpret_cast<fftwf_complex*>(outdata + i * cblock_stride);
            fftwf_execute_dft_r2c(*sforward, rdata, cdata);
        }
    }
    //! \brief Backward transform, single precision.
    void backward(std::complex<float> indata[], float outdata[], std::complex<float>*) const override{
        make_plan(sbackward);
        for(int i=0; i<blocks; i++){
            fftwf_complex* cdata = const_cast<fftwf_complex*>(reinterpret_cast<fftwf_complex const*>(indata + i * cblock_stride));
            fftwf_execute_dft_c2r(*sbackward, cdata, outdata + i * rblock_stride);
        }
    }
    //! \brief Forward transform, double precision.
    void forward(double const indata[], std::complex<double> outdata[], std::complex<double>*) const override{
        make_plan(dforward);
        for(int i=0; i<blocks; i++){
            double *rdata = const_cast<double*>(indata + i * rblock_stride);
            fftw_complex* cdata = reinterpret_cast<fftw_complex*>(outdata + i * cblock_stride);
            fftw_execute_dft_r2c(*dforward, rdata, cdata);
        }
    }
    //! \brief Backward transform, double precision.
    void backward(std::complex<double> indata[], double outdata[], std::complex<double>*) const override{
        make_plan(dbackward);
        for(int i=0; i<blocks; i++){
            fftw_complex* cdata = const_cast<fftw_complex*>(reinterpret_cast<fftw_complex const*>(indata + i * cblock_stride));
            fftw_execute_dft_c2r(*dbackward, cdata, outdata + i * rblock_stride);
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
    void make_plan(std::unique_ptr<plan_fftw<scalar_type, dir>> &plan) const{
        if (!plan) plan = std::unique_ptr<plan_fftw<scalar_type, dir>>(new plan_fftw<scalar_type, dir>(size, howmanyffts, stride, rdist, cdist));
    }

    int size, howmanyffts, stride, blocks;
    int rdist, cdist, rblock_stride, cblock_stride, rsize, csize;
    mutable std::unique_ptr<plan_fftw<float, direction::forward>> sforward;
    mutable std::unique_ptr<plan_fftw<double, direction::forward>> dforward;
    mutable std::unique_ptr<plan_fftw<float, direction::backward>> sbackward;
    mutable std::unique_ptr<plan_fftw<double, direction::backward>> dbackward;
};

/*!
 * \internal
 * \ingroup hefftefftw
 * \brief Wrapper for the single/double precision plan, make and destroy (double precision).
 *
 * \endinternal
 */
template<typename>
struct fftw_r2r_types {
    //! \brief Plan type.
    using plan_type = fftw_plan;
    //! \brief Make the plan.
    template<typename... vars> static plan_type plan_many(vars... args){ return fftw_plan_many_r2r(args...); }
    //! \brief Destructor.
    static void plan_destroy(plan_type p){ fftw_destroy_plan(p); }
};
/*!
 * \internal
 * \ingroup hefftefftw
 * \brief Wrapper for the single/double precision plan, make and destroy (single precision).
 *
 * \endinternal
 */
template<> struct fftw_r2r_types<float> {
    //! \brief Plan type.
    using plan_type = fftwf_plan;
    //! \brief Make the plan.
    template<typename... vars> static plan_type plan_many(vars... args){ return fftwf_plan_many_r2r(args...); }
    //! \brief Destructor.
    static void plan_destroy(plan_type p){ fftwf_destroy_plan(p); }
};
/*!
 * \ingroup hefftefftw
 * \brief Wrapper for the r2r plan to allow for RAII style of management.
 *
 * The plan handles scalar_type as float/double, the preprocessor defines whether to use sin/cos transform,
 * and direction.
 */
template<typename scalar_type, typename preprocessor, direction dir>
struct plan_fftw_r2r{
    //! \brief Identical to the float-complex specialization.
    plan_fftw_r2r(int size, int howmanyffts, int stride, int dist){
        auto kind = make_kind_array<1>();
        plan = fftw_r2r_types<scalar_type>::plan_many(1, &size, howmanyffts, nullptr, nullptr, stride, dist,
                                                      nullptr, nullptr, stride, dist, kind.data(), FFTW_ESTIMATE);
    }
    //! \brief Identical to the float-complex specialization.
    plan_fftw_r2r(int size1, int size2, std::array<int, 2> const &embed, int howmanyffts, int stride, int dist){
        std::array<int, 2> size = {size2, size1};
        auto kind = make_kind_array<2>();

        if (embed[0] == 0 and embed[1] == 0){
            plan = fftw_r2r_types<scalar_type>::plan_many(2, size.data(), howmanyffts, nullptr, nullptr, stride, dist,
                                                          nullptr, nullptr, stride, dist, kind.data(), FFTW_ESTIMATE);
        }else{
            plan = fftw_r2r_types<scalar_type>::plan_many(2, size.data(), howmanyffts, nullptr, embed.data(), stride, dist,
                                                          nullptr, embed.data(), stride, dist, kind.data(), FFTW_ESTIMATE);
        }
    }
    //! \brief Identical to the float-complex specialization.
    plan_fftw_r2r(int size1, int size2, int size3){
        std::array<int, 3> size = {size3, size2, size1};
        auto kind = make_kind_array<3>();
        plan = fftw_r2r_types<scalar_type>::plan_many(3, size.data(), 1, nullptr, nullptr, 1, 1, nullptr, nullptr, 1, 1, kind.data(), FFTW_ESTIMATE);
    }
    //! \brief Make the array with the kind of the transform.
    template<size_t dims> static std::array<fftw_r2r_kind, dims> make_kind_array() {
        std::array<fftw_r2r_kind, dims> kind;
        if (std::is_same<preprocessor, cpu_cos_pre_pos_processor>::value) {
            kind[0] = (dir == direction::forward) ? FFTW_REDFT10 : FFTW_REDFT01;
        } else if (std::is_same<preprocessor, cpu_sin_pre_pos_processor>::value) { // sin transform
            kind[0] = (dir == direction::forward) ? FFTW_RODFT10 : FFTW_RODFT01;
        } else if (std::is_same<preprocessor, cpu_cos1_pre_pos_processor>::value) {
            kind[0] = FFTW_REDFT00;
        } else if (std::is_same<preprocessor, cpu_sin1_pre_pos_processor>::value) { // sin transform
            kind[0] = FFTW_RODFT00;
        }
        for(size_t i=1; i<kind.size(); i++) kind[i] = kind[0];
        return kind;
    }
    //! \brief Identical to the other specialization.
    ~plan_fftw_r2r(){ fftw_r2r_types<scalar_type>::plan_destroy(plan); }
    //! \brief Automatically converts to the correct plan type.
    operator typename fftw_r2r_types<scalar_type>::plan_type() const{ return plan; }
    //! \brief The actual fftw plan.
    typename fftw_r2r_types<scalar_type>::plan_type plan;
};
template<typename prepost_processor> struct real2real_executor<backend::fftw, prepost_processor> : public executor_base{
    //! \brief Construct a plan for batch 1D transforms.
    template<typename index>
    real2real_executor(void*, box3d<index> const box, int dimension) :
        size(box.size[dimension]), size2(0),
        howmanyffts(fft1d_get_howmany(box, dimension)),
        stride(fft1d_get_stride(box, dimension)),
        dist((dimension == box.order[0]) ? size : 1),
        blocks((dimension == box.order[1]) ? box.osize(2) : 1),
        block_stride(box.osize(0) * box.osize(1)),
        total_size(box.count()),
        embed({0, 0})
    {}
    //! \brief Construct a plan for batch 2D transforms, not implemented currently.
    template<typename index>
    real2real_executor(void*, box3d<index> const box, int dir1, int dir2) :
        size(box.size[std::min(dir1, dir2)]), size2(box.size[std::max(dir1, dir2)]),
        blocks(1), block_stride(0), total_size(box.count()), embed({0, 0})
    {
        int odir1 = box.find_order(dir1);
        int odir2 = box.find_order(dir2);

        if (std::min(odir1, odir2) == 0 and std::max(odir1, odir2) == 1){
            stride = 1;
            dist = size * size2;
            howmanyffts = box.size[2];
        }else if (std::min(odir1, odir2) == 1 and std::max(odir1, odir2) == 2){
            stride = box.size[0];
            dist = 1;
            howmanyffts = box.size[0];
        }else{ // case of directions (0, 2)
            stride = 1;
            dist = size;
            embed = {static_cast<int>(box.size[2]), static_cast<int>(box.size[1] * box.size[0])};
            howmanyffts = box.size[1];
        }
    }
    //! \brief Construct a plan for a single 3D transform, not implemented currently.
    template<typename index>
    real2real_executor(void*, box3d<index> const box) :
        size(box.size[0]), size2(box.size[1]), howmanyffts(box.size[2]),
        stride(0), dist(0),
        blocks(1), block_stride(0),
        total_size(box.count()),
        embed({0, 0})
    {}

    //! \brief Returns the size of the box.
    int box_size() const override{ return total_size; }
    //! \brief Returns the size of the box.
    size_t workspace_size() const override{
        return 0;
    }
    //! \brief Forward r2r, single precision.
    virtual void forward(float data[], float*) const override{
        make_plan(sforward);
        for(int i=0; i<blocks; i++){
            fftwf_execute_r2r(*sforward, data + i * block_stride, data + i * block_stride);
        }
    }
    //! \brief Forward r2r, double precision.
    virtual void forward(double data[], double*) const override{
        make_plan(dforward);
        for(int i=0; i<blocks; i++){
            fftw_execute_r2r(*dforward, data + i * block_stride, data + i * block_stride);
        }
    }
    //! \brief Backward r2r, single precision.
    virtual void backward(float data[], float*) const override{
        make_plan(sbackward);
        for(int i=0; i<blocks; i++){
            fftwf_execute_r2r(*sbackward, data + i * block_stride, data + i * block_stride);
        }
    }
    //! \brief Backward r2r, double precision.
    virtual void backward(double data[], double*) const override{
        make_plan(dbackward);
        for(int i=0; i<blocks; i++){
            fftw_execute_r2r(*dbackward, data + i * block_stride, data + i * block_stride);
        }
    }

private:
    //! \brief Helper template to initialize the plan.
    template<typename scalar_type, direction dir>
    void make_plan(std::unique_ptr<plan_fftw_r2r<scalar_type, prepost_processor, dir>> &plan) const{
        if (not plan){
            if (dist == 0)
                plan = std::unique_ptr<plan_fftw_r2r<scalar_type, prepost_processor, dir>>(new plan_fftw_r2r<scalar_type, prepost_processor, dir>(size, size2, howmanyffts));
            else if (size2 == 0)
                plan = std::unique_ptr<plan_fftw_r2r<scalar_type, prepost_processor, dir>>(new plan_fftw_r2r<scalar_type, prepost_processor, dir>(size, howmanyffts, stride, dist));
            else
                plan = std::unique_ptr<plan_fftw_r2r<scalar_type, prepost_processor, dir>>(new plan_fftw_r2r<scalar_type, prepost_processor, dir>(size, size2, embed, howmanyffts, stride, dist));
        }
    }

    int size, size2, howmanyffts, stride, dist, blocks, block_stride, total_size;
    std::array<int, 2> embed;
    mutable std::unique_ptr<plan_fftw_r2r<float, prepost_processor, direction::forward>> sforward;
    mutable std::unique_ptr<plan_fftw_r2r<double, prepost_processor, direction::forward>> dforward;
    mutable std::unique_ptr<plan_fftw_r2r<float, prepost_processor, direction::backward>> sbackward;
    mutable std::unique_ptr<plan_fftw_r2r<double, prepost_processor, direction::backward>> dbackward;
};

/*!
 * \ingroup hefftefftw
 * \brief Helper struct that defines the types and creates instances of one-dimensional executors.
 *
 * The struct is specialized for each backend.
 */
template<> struct one_dim_backend<backend::fftw>{
    //! \brief Defines the complex-to-complex executor.
    using executor = fftw_executor;
    //! \brief Defines the real-to-complex executor.
    using executor_r2c = fftw_executor_r2c;
};

/*!
 * \ingroup hefftefftw
 * \brief Helper struct that defines the types and creates instances of one-dimensional executors.
 *
 * The struct is specialized for each backend.
 */
template<> struct one_dim_backend<backend::fftw_cos>{
    //! \brief Defines the real-to-real executor.
    using executor = real2real_executor<backend::fftw, cpu_cos_pre_pos_processor>;
    //! \brief There is no real-to-complex variant.
    using executor_r2c = void;
};
/*!
 * \ingroup hefftefftw
 * \brief Helper struct that defines the types and creates instances of one-dimensional executors.
 *
 * The struct is specialized for each backend.
 */
template<> struct one_dim_backend<backend::fftw_sin>{
    //! \brief Defines the real-to-real executor.
    using executor = real2real_executor<backend::fftw, cpu_sin_pre_pos_processor>;
    //! \brief There is no real-to-complex variant.
    using executor_r2c = void;
};
/*!
 * \ingroup hefftefftw
 * \brief Helper struct that defines the types and creates instances of one-dimensional executors.
 *
 * The struct is specialized for each backend.
 */
template<> struct one_dim_backend<backend::fftw_cos1>{
    //! \brief Defines the real-to-real executor.
    using executor = real2real_executor<backend::fftw, cpu_cos1_pre_pos_processor>;
    //! \brief There is no real-to-complex variant.
    using executor_r2c = void;
};
/*!
 * \ingroup hefftefftw
 * \brief Helper struct that defines the types and creates instances of one-dimensional executors.
 *
 * The struct is specialized for each backend.
 */
template<> struct one_dim_backend<backend::fftw_sin1>{
    //! \brief Defines the real-to-real executor.
    using executor = real2real_executor<backend::fftw, cpu_sin1_pre_pos_processor>;
    //! \brief There is no real-to-complex variant.
    using executor_r2c = void;
};

/*!
 * \ingroup hefftefftw
 * \brief Sets the default options for the fftw backend.
 */
template<> struct default_plan_options<backend::fftw>{
    //! \brief The reshape operations will also reorder the data.
    static const bool use_reorder = true;
};

/*!
 * \ingroup hefftefftw
 * \brief Sets the default options for the fftw backend.
 */
template<> struct default_plan_options<backend::fftw_cos>{
    //! \brief The reshape operations will also reorder the data.
    static const bool use_reorder = true;
};
/*!
 * \ingroup hefftefftw
 * \brief Sets the default options for the fftw backend.
 */
template<> struct default_plan_options<backend::fftw_sin>{
    //! \brief The reshape operations will also reorder the data.
    static const bool use_reorder = true;
};
/*!
 * \ingroup hefftefftw
 * \brief Sets the default options for the fftw backend.
 */
template<> struct default_plan_options<backend::fftw_cos1>{
    //! \brief The reshape operations will also reorder the data.
    static const bool use_reorder = true;
};
/*!
 * \ingroup hefftefftw
 * \brief Sets the default options for the fftw backend.
 */
template<> struct default_plan_options<backend::fftw_sin1>{
    //! \brief The reshape operations will also reorder the data.
    static const bool use_reorder = true;
};

}

#endif

#endif   /* HEFFTE_BACKEND_FFTW_H */
