/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_BACKEND_STOCK_FFT_H
#define HEFFTE_BACKEND_STOCK_FFT_H

#include "heffte_r2r_executor.h"

#include "stock_fft/heffte_stock_tree.h"

/*!
 * \ingroup fft3d
 * \addtogroup hefftestock Backend stock
 *
 * Wrappers and template specializations related to the stock backend.
 */

namespace heffte{

namespace backend{
    /*!
     * \ingroup hefftestock
     * \brief Indicate that the stock backend has been enabled.
     */
    template<> struct is_enabled<stock> : std::true_type{};
    /*!
     * \ingroup hefftestock
     * \brief Indicate that the stock backend has been enabled.
     */
    template<> struct is_enabled<stock_cos> : std::true_type{};
    /*!
     * \ingroup hefftestock
     * \brief Indicate that the stock backend has been enabled.
     */
    template<> struct is_enabled<stock_sin> : std::true_type{};

// Specialization is not necessary since the default behavior assumes CPU parameters.
//     template<>
//     struct buffer_traits<stock>{
//         using location = tag::cpu;
//         template<typename T> using container = std::vector<T>;
//     };
}

//! \brief Recognize stock FFT single complex (which are std::complex) types
template<> struct is_ccomplex<stock::Complex<float, 1>> : std::true_type{};
template<> struct is_zcomplex<stock::Complex<double, 1>> : std::true_type{};
#ifdef Heffte_ENABLE_AVX
/*!
 * \ingroup hefftestock
 * \brief Recognize the stock FFT single precision complex types.
 */
template<> struct is_ccomplex<stock::Complex<float, 4>> : std::true_type{};
template<> struct is_ccomplex<stock::Complex<float, 8>> : std::true_type{};
#ifdef Heffte_ENABLE_AVX512
template<> struct is_ccomplex<stock::Complex<float, 16>> : std::true_type{};
#endif // Heffte_ENABLE_AVX512

/*!
 * \ingroup hefftestock
 * \brief Recognize the stock FFT double precision complex types.
 */
template<> struct is_zcomplex<stock::Complex<double, 2>> : std::true_type{};
template<> struct is_zcomplex<stock::Complex<double, 4>> : std::true_type{};
#ifdef Heffte_ENABLE_AVX512
template<> struct is_zcomplex<stock::Complex<double, 8>> : std::true_type{};
#endif // Heffte_ENABLE_AVX512
#endif // Heffte_ENABLE_AVX

#ifdef Heffte_ENABLE_AVX

//! \brief Copy an array of numbers into a stock::Complex where only the first c_len spots are filled
template<typename F, int L>
stock::Complex<F, L> copy_pad(std::complex<F> const *c, int c_len, int i_stride) {
    if(c_len > L/2) {
        throw std::runtime_error("Invalid padding requested!\n");
    }
    std::complex<F> ret[L];
    for(int i = 0; i < c_len; i++) ret[i] = c[i*i_stride];
    return stock::Complex<F,L> (ret);
}
//! \brief Copy an array of numbers into a stock::Complex where only the first c_len spots are filled
template<typename F, int L>
stock::Complex<F,L> copy_pad(F const *c, int c_len, int i_stride) {
    if(c_len > L/2) {
        throw std::runtime_error("Invalid padding requested!\n");
    }
    std::complex<F> ret[L];
    for(int i = 0; i < c_len; i++) ret[i] = std::complex<F>(c[i*i_stride], 0.0);
    return stock::Complex<F,L> (ret);
}

template<typename F> struct pack_size { };
#ifdef Heffte_ENABLE_AVX512
template<>           struct pack_size<float> {constexpr static int size = 16;};
template<>           struct pack_size<double>{constexpr static int size = 8;};
#else
template<>           struct pack_size<float> {constexpr static int size = 8;};
template<>           struct pack_size<double>{constexpr static int size = 4;};
#endif // Heffte_ENABLE_AVX512

/*!
 * \ingroup hefftestock
 * \brief Base plan for stock library, using only the specialization for float and double complex.
 *
 * Stock fft uses plans for forward and backward fft transforms. The default implementation is
 * for r2c transformations
 * The specializations to this struct will wrap around such plans and provide RAII style
 * of memory management and simple constructors that take inputs suitable to HeFFTe.
 */
template<typename F, direction dir>
struct plan_stock_fft{
    /*!
     * \brief Constructor taking into account the different sizes for the real and complex parts.
     *
     * \param size is the number of entries in a 1-D transform
     * \param howmanyffts is the number of transforms in the batch
     * \param stride is the distance between entries of the same transform
     * \param rdist is the distance between the first entries of consecutive sequences in the real sequences
     * \param cdist is the distance between the first entries of consecutive sequences in the complex sequences
     */
    plan_stock_fft(int size, int howmanyffts, int stride, int rdist, int cdist):
                    N(size), num_ffts(howmanyffts), stride_sz(stride), real_d(rdist), comp_d(cdist) {
        numNodes = stock::getNumNodes(N);
        root = std::unique_ptr<stock::biFuncNode<F,L>[]>(new stock::biFuncNode<F,L>[numNodes]);
        init_fft_tree(root.get(), N);
    }
    //! \brief Identical to the F-complex specialization.
    int N, num_ffts, stride_sz, real_d, comp_d, numNodes;
    constexpr static int L = pack_size<F>::size;
    std::unique_ptr<stock::biFuncNode<F,L>[]> root;
    //! \brief Execute C2R FFT
    void execute(std::complex<F> const idata[], F odata[]) {
        // Allocate input and output temporary arrays
        stock::complex_vector<F,L> inp (N);
        stock::complex_vector<F,L> out (N);

        // Perform batch transform on everything save for the remainder
        for(int p = 0; p+((L/2)-1) < num_ffts; p += (L/2)) {
            // Convert types
            inp[0] = stock::Complex<F,L> (&idata[p*comp_d], comp_d);
            for(int i = 1; i < (N+2)/2; i++) {
                int idx = p*comp_d + i*stride_sz;
                inp[i] = stock::Complex<F,L> (&idata[idx], comp_d);
                inp[N-i] = inp[i].conjugate();
            }
            // Perform fft
            root[0].fptr(inp.data(), out.data(), 1, 1, root.get(), dir);
            // Convert type back
            for(int i = 0; i < N; i++) {
                int idx = p*real_d + i*stride_sz;
                stock::Complex<F,L> arr = out[i];
                for(int j = 0; j < L/2; j++) odata[idx+j*real_d] = arr[j].real();
            }
        }

        // Handle remainder
        if(num_ffts % (L/2) > 0) {
            int rem = num_ffts % (L/2);
            // Init p for ease of use
            int p = num_ffts - rem;
            inp[0] = copy_pad<F,L>(&idata[p*comp_d], rem, comp_d);
            for(int i = 1; i < (N+2)/2; i++) {
                int idx = p*comp_d + i*stride_sz;
                inp[i] = copy_pad<F,L>(&idata[idx], rem, comp_d);
                inp[N-i] = inp[i].conjugate();
            }
            root[0].fptr(inp.data(), out.data(), 1, 1, root.get(), dir);
            for(int i = 0; i < N; i++) {
                int idx = p*real_d + i*stride_sz;
                stock::Complex<F,L> arr = out[i];
                // Only need first k columns
                for(int k = 0; k < rem; k++) {
                    odata[idx + k*real_d] = arr[k].real();
                }
            }
        }
    }
    //! \brief Execute R2C FFT
    void execute(F const idata[], std::complex<F> odata[]) {
        // Allocate input and output temporary arrays
        stock::complex_vector<F,L> inp (N);
        stock::complex_vector<F,L> out (N);

        // Perform batch transform on everything save for the remainder
        for(int p = 0; p+((L/2)-1) < num_ffts; p += L/2) {
            // Convert types
            for(int i = 0; i < N; i++) {
                int idx = p*real_d + i*stride_sz;
                inp[i] = copy_pad<F,L>(&idata[idx], L/2, real_d);
            }
            // Perform fft
            root[0].fptr(inp.data(), out.data(), 1, 1, root.get(), dir);
            // Convert type back
            for(int i = 0; i < (N+2)/2; i++) {
                int idx = p*comp_d + i*stride_sz;
                stock::Complex<F,L> arr = out[i];
                for(int j = 0; j < L/2; j++) odata[idx+j*comp_d] = arr[j];
            }
        }

        // Handle remainder
        if(num_ffts % (L/2) > 0) {
            int rem = num_ffts % (L/2);
            // Init p for ease of use
            int p = num_ffts - rem;
            for(int i = 0; i < N; i++) {
                int idx = p*real_d + i*stride_sz;
                // remainder columns are all zeros
                inp[i] = copy_pad<F,L>(&idata[idx], rem, real_d);
            }
            root[0].fptr(inp.data(), out.data(), 1, 1, root.get(), dir);
            for(int i = 0; i < (N+2)/2; i++) {
                int idx = p*comp_d + i*stride_sz;
                stock::Complex<F,L> arr = out[i];
                // Only need first k columns
                for(int k = 0; k < rem; k++) {
                    odata[idx + k*comp_d] = arr[k];
                }
            }
        }
    }
};

/*!
 * \ingroup hefftestock
 * \brief Plan for the single or double precision complex transform.
 *
 * \param dir indicates a forward or backward transform
 */
template<typename F, direction dir>
struct plan_stock_fft<std::complex<F>, dir>{
    /*!
     * \brief Constructor to plan out an FFT using the stock backend
     *
     * \param size is the number of entries in a 1-D transform
     * \param howmanyffts is the number of transforms in the batch
     * \param stride is the distance between entries of the same transform
     * \param dist is the distance between the first entries of consecutive sequences
     */
    plan_stock_fft(int size, int howmanyffts, int stride, int dist):
                   N(size), num_ffts(howmanyffts), stride_sz(stride), idist(dist), odist(dist) {
        numNodes = stock::getNumNodes(N);
        root = std::unique_ptr<stock::biFuncNode<F,L>[]>(new stock::biFuncNode<F,L>[numNodes]);
        init_fft_tree(root.get(), N);
    }

    int N, num_ffts, stride_sz, idist, odist, numNodes;
    constexpr static int L = pack_size<F>::size;
    std::unique_ptr<stock::biFuncNode<F, L>[]> root;
    //! \brief Execute an FFT inplace on std::complex<F> data
    void execute(std::complex<F> data[]) {
        // Allocate input and output temporary arrays
        stock::complex_vector<F,L> inp (N);
        stock::complex_vector<F,L> out (N);

        // Perform batch transform on everything save for the remainder
        for(int p = 0; p + (L/2 - 1) < num_ffts; p += L/2) {
            // Convert types
            for(int i = 0; i < N; i++) {
                int idx = p*idist + i*stride_sz;
                inp[i] = stock::Complex<F,L> (&data[idx], idist);
            }
            // Perform fft
            root[0].fptr(inp.data(), out.data(), 1, 1, root.get(), dir);
            // Convert type back
            for(int i = 0; i < N; i++) {
                int idx = p*odist + i*stride_sz;
                stock::Complex<F,L> arr = out[i];
                for(int j = 0; j < L/2; j++) data[idx+j*odist] = arr[j];
            }
        }

        // Handle remainder
        if(num_ffts % (L/2) > 0) {
            int rem = num_ffts % (L/2);
            // Init p for ease of use
            int p = num_ffts - rem;
            for(int i = 0; i < N; i++) {
                int idx = p*idist + i*stride_sz;
                // remainder columns are all zeros
                inp[i] = copy_pad<F,L>(&data[idx], rem, idist);
            }
            root[0].fptr(inp.data(), out.data(), 1, 1, root.get(), dir);
            for(int i = 0; i < N; i++) {
                int idx = p*odist + i*stride_sz;
                stock::Complex<F,L> arr = out[i];
                // Only need first k columns
                for(int k = 0; k < rem; k++) {
                    data[idx + k*odist] = arr[k];
                }
            }
        }
    }
};

#else

// In the case the user does not have AVX installed

/*!
 * \ingroup hefftestock
 * \brief Specialization for r2c single precision.
 */
template<typename F, direction dir>
struct plan_stock_fft{
    /*!
     * \brief Constructor taking into account the different sizes for the real and complex parts.
     *
     * \param size is the number of entries in a 1-D transform
     * \param howmanyffts is the number of transforms in the batch
     * \param stride is the distance between entries of the same transform
     * \param rdist is the distance between the first entries of consecutive sequences in the real sequences
     * \param cdist is the distance between the first entries of consecutive sequences in the complex sequences
     */
    plan_stock_fft(int size, int howmanyffts, int stride, int rdist, int cdist):
                    N(size), num_ffts(howmanyffts), stride_sz(stride), real_d(rdist), comp_d(cdist) {
        numNodes = stock::getNumNodes(N);
        root = std::unique_ptr<stock::biFuncNode<F,1>[]>(new stock::biFuncNode<F,1>[numNodes]);
        init_fft_tree(root.get(), N);
    }
    //! \brief Identical to the F-complex specialization.
    int N, num_ffts, stride_sz, real_d, comp_d, numNodes;
    std::unique_ptr<stock::biFuncNode<F, 1>[]> root;
    //! \brief Execute C2R FFT
    void execute(std::complex<F> const idata[], F odata[]) {
        // Allocate input and output temporary arrays
        stock::complex_vector<F,1> inp (N);
        stock::complex_vector<F,1> out (N);

        // Perform batch transform on everything save for the remainder
        for(int p = 0; p < num_ffts; p++) {
            // Convert types
            inp[0] = stock::Complex<F,1> {idata[p*comp_d]};
            for(int i = 1; i < (N+2)/2; i++) {
                int idx = p*comp_d + i*stride_sz;
                inp[i] = stock::Complex<F,1> {idata[idx]};
                inp[N-i] = inp[i].conjugate();
            }
            // Perform fft
            root[0].fptr(inp.data(), out.data(), 1, 1, root.get(), dir);
            // Convert type back
            for(int i = 0; i < N; i++) {
                int idx = p*real_d + i*stride_sz;
                stock::Complex<F,1> arr = out[i];
                odata[idx+0*real_d] = arr[0].real();
            }
        }
    }
    //! \brief Execute R2C FFT
    void execute(F const idata[], std::complex<F> odata[]) {
        // Allocate input and output temporary arrays
        stock::complex_vector<F,1> inp (N);
        stock::complex_vector<F,1> out (N);

        // Perform batch transform on everything save for the remainder
        for(int p = 0; p < num_ffts; p++) {
            // Convert types
            for(int i = 0; i < N; i++) {
                int idx = p*real_d + i*stride_sz;
                inp[i] = stock::Complex<F,1> {idata[idx+0*real_d], 0};
            }
            // Perform fft
            root[0].fptr(inp.data(), out.data(), 1, 1, root.get(), dir);
            // Convert type back
            for(int i = 0; i < (N+2)/2; i++) {
                int idx = p*comp_d + i*stride_sz;
                stock::Complex<F,1> arr = out[i];
                odata[idx+0*comp_d] = arr[0];
            }
        }
    }
};

/*!
 * \ingroup hefftestock
 * \brief Plan for the single precision complex transform.
 *
 * \tparam dir indicates a forward or backward transform
 */
template<typename F, direction dir>
struct plan_stock_fft<std::complex<F>, dir>{
    /*!
     * \brief Constructor to plan out an FFT using the stock backend
     *
     * \param size is the number of entries in a 1-D transform
     * \param howmanyffts is the number of transforms in the batch
     * \param stride is the distance between entries of the same transform
     * \param dist is the distance between the first entries of consecutive sequences
     */
    plan_stock_fft(int size, int howmanyffts, int stride, int dist):
                    N(size), num_ffts(howmanyffts), stride_sz(stride), idist(dist), odist(dist) {
        numNodes = stock::getNumNodes(N);
        root = std::unique_ptr<stock::biFuncNode<F,1>[]>(new stock::biFuncNode<F,1>[numNodes]);
        init_fft_tree(root.get(), N);
    }

    int N, num_ffts, stride_sz, idist, odist, numNodes;
    std::unique_ptr<stock::biFuncNode<F,1>[]> root;
    //! \brief Execute an FFT inplace on std::complex<F> data
    void execute(std::complex<F> data[]) {
        // Allocate input and output temporary arrays
        stock::complex_vector<F,1> inp (N);
        stock::complex_vector<F,1> out (N);

        // Perform batch transform on everything save for the remainder
        for(int p = 0; p < num_ffts; p++) {
            // Convert types
            for(int i = 0; i < N; i++) {
                int idx = p*idist + i*stride_sz;
                inp[i] = stock::Complex<F, 1> {data[idx+0*idist]};
            }
            // Perform fft
            root[0].fptr(inp.data(), out.data(), 1, 1, root.get(), dir);
            // Convert type back
            for(int i = 0; i < N; i++) {
                int idx = p*odist + i*stride_sz;
                stock::Complex<F, 1> arr = out[i];
                data[idx+0*odist] = arr[0];
            }
        }
    }
};
#endif // Heffte_ENABLE_AVX512

/*!
 * \ingroup hefftestock
 * \brief Wrapper around the Stock FFT API.
 *
 * A single class that manages the plans and executions of stock fft
 * so that a single API is provided for all backends.
 * The executor operates on a box and performs 1-D FFTs
 * for the given dimension.
 * The class silently manages the plans and buffers needed
 * for the different types.
 * All input and output arrays must have size equal to the box.
 */
class stock_fft_executor : public executor_base{
public:
    //! \brief Bring forth method that have not been overloaded.
    using executor_base::forward;
    //! \brief Bring forth method that have not been overloaded.
    using executor_base::backward;
    //! \brief Bring forth method that have not been overloaded.
    using executor_base::complex_size;
    //! \brief Constructor, specifies the box and dimension.
    template<typename index>
    stock_fft_executor(void*, box3d<index> const box, int dimension) :
        size(box.size[dimension]),
        num_ffts(fft1d_get_howmany(box, dimension)),
        stride(fft1d_get_stride(box, dimension)),
        dist((dimension == box.order[0]) ? size : 1),
        blocks((dimension == box.order[1]) ? box.osize(2) : 1),
        block_stride(box.osize(0) * box.osize(1)),
        total_size(box.count())
    {}
    //! \brief Placeholder, unimplemented.
    template<typename index> stock_fft_executor(void*, box3d<index> const, int, int){
        throw std::runtime_error("2D transforms for the stock backend are not available yet!");
    }
    //! \brief Placeholder, unimplemented.
    template<typename index> stock_fft_executor(void*, box3d<index> const){
        throw std::runtime_error("3D transforms for the stock backend are not available yet!");
    }

    //! \brief Forward fft, float-complex case.
    void forward(std::complex<float> data[], std::complex<float>*) const override{
        make_plan(cforward);
        for(int i=0; i<blocks; i++) {
            cforward->execute(data + i * block_stride);
        }
    }
    //! \brief Backward fft, float-complex case.
    void backward(std::complex<float> data[], std::complex<float>*) const override{
        make_plan(cbackward);
        for(int i=0; i<blocks; i++) {
            cbackward->execute(data + i * block_stride);
        }
    }
    //! \brief Forward fft, double-complex case.
    void forward(std::complex<double> data[], std::complex<double>*) const override{
        make_plan(zforward);
        for(int i=0; i<blocks; i++) {
            zforward->execute(data + i * block_stride);
        }
    }
    //! \brief Backward fft, double-complex case.
    void backward(std::complex<double> data[], std::complex<double>*) const override{
        make_plan(zbackward);
        for(int i=0; i<blocks; i++) {
            zbackward->execute(data + i * block_stride);
        }
    }

    //! \brief Converts the real data to complex and performs float-complex forward transform.
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
    void make_plan(std::unique_ptr<plan_stock_fft<scalar_type, dir>> &plan) const{
        if (!plan) plan = std::unique_ptr<plan_stock_fft<scalar_type, dir>>(new plan_stock_fft<scalar_type, dir>(size, num_ffts, stride, dist));
    }

    int size, num_ffts, stride, dist, blocks, block_stride, total_size;
    mutable std::unique_ptr<plan_stock_fft<std::complex<float>, direction::forward>> cforward;
    mutable std::unique_ptr<plan_stock_fft<std::complex<float>, direction::backward>> cbackward;
    mutable std::unique_ptr<plan_stock_fft<std::complex<double>, direction::forward>> zforward;
    mutable std::unique_ptr<plan_stock_fft<std::complex<double>, direction::backward>> zbackward;
};

/*!
 * \ingroup hefftestock
 * \brief Wrapper to stock API for real-to-complex transform with shortening of the data.
 *
 * Serves the same purpose of heffte::stock_fft_executor but only real input is accepted
 * and only the unique (non-conjugate) coefficients are computed.
 * All real arrays must have size of real_size() and all complex arrays must have size complex_size().
 */
class stock_fft_executor_r2c : public executor_base{
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
    stock_fft_executor_r2c(void*, box3d<index> const box, int dimension) :
        size(box.size[dimension]),
        num_ffts(fft1d_get_howmany(box, dimension)),
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
        for(int i=0; i<blocks; i++) {
            sforward->execute(indata + i*rblock_stride, outdata + i*cblock_stride);
        }
    }
    //! \brief Backward transform, single precision.
    void backward(std::complex<float> indata[], float outdata[], std::complex<float>*) const override{
        make_plan(sbackward);
        for(int i=0; i<blocks; i++) {
            sbackward->execute(indata + i*cblock_stride, outdata + i*rblock_stride);
        }
    }
    //! \brief Forward transform, double precision.
    void forward(double const indata[], std::complex<double> outdata[], std::complex<double>*) const override{
        make_plan(dforward);
        for(int i=0; i<blocks; i++) {
            dforward->execute(indata + i*rblock_stride, outdata + i*cblock_stride);
        }
    }
    //! \brief Backward transform, double precision.
    void backward(std::complex<double> indata[], double outdata[], std::complex<double>*) const override{
        make_plan(dbackward);
        for(int i=0; i<blocks; i++) {
            dbackward->execute(indata + i*cblock_stride, outdata + i*rblock_stride);
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
    void make_plan(std::unique_ptr<plan_stock_fft<scalar_type, dir>> &plan) const{
        if (!plan) plan = std::unique_ptr<plan_stock_fft<scalar_type, dir>>(new plan_stock_fft<scalar_type, dir>(size, num_ffts, stride, rdist, cdist));
    }

    int size, num_ffts, stride, blocks;
    int rdist, cdist, rblock_stride, cblock_stride, rsize, csize;
    mutable std::unique_ptr<plan_stock_fft<float, direction::forward>> sforward;
    mutable std::unique_ptr<plan_stock_fft<double, direction::forward>> dforward;
    mutable std::unique_ptr<plan_stock_fft<float, direction::backward>> sbackward;
    mutable std::unique_ptr<plan_stock_fft<double, direction::backward>> dbackward;
};
/*!
 * \ingroup hefftestock
 * \brief Helper struct that defines the types and creates instances of one-dimensional executors.
 *
 * The struct is specialized for each backend.
 */
template<> struct one_dim_backend<backend::stock>{
    //! \brief Defines the complex-to-complex executor.
    using executor = stock_fft_executor;
    //! \brief Defines the real-to-complex executor.
    using executor_r2c = stock_fft_executor_r2c;
};
/*!
 * \ingroup hefftestock
 * \brief Helper struct that defines the types and creates instances of one-dimensional executors.
 *
 * The struct is specialized for each backend.
 */
template<> struct one_dim_backend<backend::stock_cos>{
    //! \brief Defines the real-to-real executor.
    using executor = real2real_executor<backend::stock, cpu_cos_pre_pos_processor>;
    //! \brief There is no real-to-complex variant.
    using executor_r2c = void;
};
/*!
 * \ingroup hefftestock
 * \brief Helper struct that defines the types and creates instances of one-dimensional executors.
 *
 * The struct is specialized for each backend.
 */
template<> struct one_dim_backend<backend::stock_sin>{
    //! \brief Defines the real-to-real executor.
    using executor = real2real_executor<backend::stock, cpu_sin_pre_pos_processor>;
    //! \brief There is no real-to-complex variant.
    using executor_r2c = void;
};

/*!
 * \ingroup hefftestock
 * \brief Sets the default options for the stock fft backend.
 */
template<> struct default_plan_options<backend::stock>{
    //! \brief The reshape operations will also reorder the data.
    static const bool use_reorder = true;
};
/*!
 * \ingroup hefftestock
 * \brief Sets the default options for the stock fft backend.
 */
template<> struct default_plan_options<backend::stock_cos>{
    //! \brief The reshape operations will also reorder the data.
    static const bool use_reorder = true;
};
/*!
 * \ingroup hefftestock
 * \brief Sets the default options for the stock fft backend.
 */
template<> struct default_plan_options<backend::stock_sin>{
    //! \brief The reshape operations will also reorder the data.
    static const bool use_reorder = true;
};

}

#endif   /* HEFFTE_BACKEND_STOCK_FFT_H */
