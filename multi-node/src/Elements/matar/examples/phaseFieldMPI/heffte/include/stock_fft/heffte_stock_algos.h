/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_STOCK_ALGOS_H
#define HEFFTE_STOCK_ALGOS_H

#include <cmath>

#include "heffte_stock_complex.h"
#include "heffte_stock_allocator.h"
#include "heffte_common.h"

//! \brief If the signal is smaller than this, we use the DFT
#define HEFFTE_STOCK_THRESHOLD 1

namespace heffte {
//! \brief Find the sign given a direction
inline int direction_sign(direction dir) {
    return (dir == direction::forward) ? -1 : 1;
}

namespace stock {
// Need forward declaration for the using directive
template<typename F, int L>
struct biFuncNode;

// Functions in the algos implementation file facing externally
/*! \ingroup stockbackend
 * \brief Create a stock Complex representation of a twiddle factor
 */
template<typename F, int L>
struct omega {
    static inline Complex<F, L> get(size_t power, size_t N, direction dir) {
        F a = 2.*M_PI*((F) power)/((F) N);
        return Complex<F,L>(cos(a), direction_sign(dir)*sin(a));
    }
};

//! \brief Enumeration to represent which FFT type to use
enum fft_type {pow2, pow3, pow4, composite, discrete, rader};

/*! \ingroup stockbackend
 * \brief Functor class to represent any Fourier Transform
 * A class to use lambdas to switch between what FFT should be called
 */
template<typename F, int L>
class Fourier_Transform {
    protected:
        fft_type type;
        size_t root = 0;
        size_t root_inv = 0;
    public:
        explicit Fourier_Transform(fft_type fft): type(fft) {}
        explicit Fourier_Transform(size_t a, size_t ainv): type(fft_type::rader), root(a), root_inv(ainv) { }
        void operator()(Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, biFuncNode<F,L>* sRoot, direction dir) {
            switch(type) {
                case fft_type::pow2: pow2_FFT(x, y, s_in, s_out, sRoot, dir); break;
                case fft_type::pow3: pow3_FFT(x, y, s_in, s_out, sRoot, dir); break;
                case fft_type::pow4: pow4_FFT(x, y, s_in, s_out, sRoot, dir); break;
                case fft_type::composite: composite_FFT(x, y, s_in, s_out, sRoot, dir); break;
                case fft_type::discrete: DFT(x, y, s_in, s_out, sRoot, dir); break;
                case fft_type::rader: rader_FFT(x, y, s_in, s_out, sRoot, dir, root, root_inv); break;
                default: throw std::runtime_error("Invalid Fourier Transform Form!\n");
            }
        }
};

/*! \ingroup stockbackend
 * \brief Class to represent the call-graph nodes
 *
 * A class to map out what the call-graph of the FFT will look like, then
 * hold it in memory for the FFTs.
 */
template<typename F, int L>
struct biFuncNode {
    size_t sz = 0;               // Size of FFT
    Fourier_Transform<F,L> fptr; // FFT for this call
    size_t left = 0;             // Offset in array until left child
    size_t right = 0;            // Offset in array until right child
    complex_vector<F,L> workspace; // Workspace
    biFuncNode(): fptr(fft_type::discrete) {};
    biFuncNode(size_t sig_size, fft_type type): sz(sig_size), fptr(type), workspace(sig_size) {}; // Create default constructor
    biFuncNode(size_t sig_size, size_t a, size_t ainv): fptr(a,ainv), workspace(sig_size) {};
};

// Internal helper function to perform a DFT
template<typename F, int L>
inline void DFT_helper(size_t size, Complex<F,L>* sig_in, Complex<F,L>* sig_out, size_t s_in, size_t s_out, direction dir) {
    if(size == 1) {
        sig_out[0] = sig_in[0];
        return;
    }

    // Twiddle with smallest numerator
    Complex<F,L> w0 = omega<F,L>::get(1, size, dir);

    // Base twiddle for each outer iteration
    Complex<F,L> wk = w0;

    // Twiddle for inner iterations
    Complex<F,L> wkn = w0;

    // Calculate first element of output
    Complex<F,L> tmp = sig_in[0];
    for(size_t n = 1; n < size; n++) {
        tmp += sig_in[n*s_in];
    }
    sig_out[0] = tmp;

    // Initialize rest of output
    for(size_t k = 1; k < size; k++) {
        // Initialize kth output
        tmp = sig_in[0];

        // Calculate kth output
        for(size_t n = 1; n < size; n++) {
            tmp = wkn.fmadd(sig_in[n*s_in], tmp);
            wkn *= wk;
        }
        sig_out[k*s_out] = tmp;

        // "Increment" wk and "reset" wkn
        wk *= w0;
        wkn = wk;
    }
}

// External-facing function to properly call the internal DFT function
template<typename F, int L>
inline void DFT(Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, biFuncNode<F,L>* sLeaf, direction dir) {
    DFT_helper(sLeaf->sz, x, y, s_in, s_out, dir);
}

// Recursive helper function implementing a classic C-T FFT
template<typename F, int L>
inline void pow2_FFT_helper(size_t N, Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, direction dir) {
    // Trivial case
    if(N == 2) {
        y[    0] = x[0] + x[s_in];
        y[s_out] = x[0] - x[s_in];
        return;
    }

    // Size of sub-problem
    int m = N/2;

    // Divide into two sub-problems
    pow2_FFT_helper(m, x, y, s_in*2, s_out, dir);
    pow2_FFT_helper(m, x+s_in, y+s_out*m, s_in*2, s_out, dir);

    // Twiddle Factor
    Complex<F,L> w1 = omega<F,L>::get(1, N, dir);
    Complex<F,L> wj = w1;
    Complex<F,L> y_j = y[0];

    y[0] += y[m*s_out];
    y[m*s_out] = y_j - y[m*s_out];

    // Conquer larger problem accordingly
    for(int j = 1; j < m; j++) {
        int j_stride = j*s_out;
        int jm_stride = (j+m)*s_out;
        y_j = y[j_stride];
        y[j_stride]  = y_j + wj*y[jm_stride];
        y[jm_stride] = y_j - wj*y[jm_stride];
        wj *= w1;
    }
}

// External function to call the C-T radix-2 FFT
template<typename F, int L>
inline void pow2_FFT(Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, biFuncNode<F,L>* sRoot, direction dir) {
    const size_t N = sRoot->sz; // Size of problem
    pow2_FFT_helper(N, x, y, s_in, s_out, dir); // Call the radix-2 FFT
}

// Recursive helper function implementing a classic C-T FFT
template<typename F, int L>
inline void pow4_FFT_helper(size_t N, Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, direction dir) {
    // Trivial case
    if(N == 4) {
        if(dir == direction::forward) {
            y[0*s_out] = x[0] + x[2*s_in] + (x[s_in] + x[3*s_in]);
            y[1*s_out] = x[0] - x[2*s_in] + (x[s_in] - x[3*s_in]).__mul_neg_i();
            y[2*s_out] = x[0] + x[2*s_in] - (x[s_in] + x[3*s_in]);
            y[3*s_out] = x[0] - x[2*s_in] + (x[s_in] - x[3*s_in]).__mul_i();
        } else {
            y[0*s_out] = x[0] + x[2*s_in] + (x[s_in] + x[3*s_in]);
            y[1*s_out] = x[0] - x[2*s_in] + (x[s_in] - x[3*s_in]).__mul_i();
            y[2*s_out] = x[0] + x[2*s_in] - (x[s_in] + x[3*s_in]);
            y[3*s_out] = x[0] - x[2*s_in] + (x[s_in] - x[3*s_in]).__mul_neg_i();
        }
        return;
    }

    // Size of sub-problem
    int m = N/4;

    // Divide into two sub-problems
    pow4_FFT_helper(m, x         , y            , s_in*4, s_out, dir);
    pow4_FFT_helper(m, x +   s_in, y +   s_out*m, s_in*4, s_out, dir);
    pow4_FFT_helper(m, x + 2*s_in, y + 2*s_out*m, s_in*4, s_out, dir);
    pow4_FFT_helper(m, x + 3*s_in, y + 3*s_out*m, s_in*4, s_out, dir);

    // Twiddle Factors
    Complex<F,L> w1 = omega<F,L>::get(1,N,dir);
    Complex<F,L> w2 = w1*w1;
    Complex<F,L> w3 = w2*w1;
    Complex<F,L> wk1 = w1; Complex<F,L> wk2 = w2; Complex<F,L> wk3 = w3;
    int k0 =         0;
    int k1 =   m*s_out;
    int k2 = 2*m*s_out;
    int k3 = 3*m*s_out;
    Complex<F,L> y_k0 = y[k0];
    Complex<F,L> y_k1 = y[k1];
    Complex<F,L> y_k2 = y[k2];
    Complex<F,L> y_k3 = y[k3];
    // Conquer larger problem accordingly
    if(dir == direction::forward) {
        y[k0] = y_k0 + y_k2 + (y_k1 + y_k3);
        y[k1] = y_k0 - y_k2 + (y_k1 - y_k3).__mul_neg_i();
        y[k2] = y_k0 + y_k2 - (y_k1 + y_k3);
        y[k3] = y_k0 - y_k2 + (y_k1 - y_k3).__mul_i();
        for(int k = 1; k < m; k++) {
            k0 = (k      )*s_out;
            k1 = (k +   m)*s_out;
            k2 = (k + 2*m)*s_out;
            k3 = (k + 3*m)*s_out;
            y_k0  = y[k0];
            y_k1  = y[k1];
            y_k2  = y[k2];
            y_k3  = y[k3];
            y[k0] = wk2.fmadd( y_k2, y_k0) + wk1.fmadd(y_k1, wk3*y_k3);
            y[k1] = wk2.fmadd(-y_k2, y_k0) + wk3.fmsub(y_k3, wk1*y_k1).__mul_i();
            y[k2] = wk2.fmadd( y_k2, y_k0) - wk1.fmadd(y_k1, wk3*y_k3);
            y[k3] = wk2.fmadd(-y_k2, y_k0) + wk1.fmsub(y_k1, wk3*y_k3).__mul_i();
            wk1 *= w1; wk2 *= w2; wk3 *= w3;
        }
    }
    else {
        y[k0] = y_k0 + y_k2 +  y_k1 + y_k3;
        y[k1] = y_k0 - y_k2 + (y_k1 - y_k3).__mul_i();
        y[k2] = y_k0 + y_k2 -  y_k1 - y_k3;
        y[k3] = y_k0 - y_k2 + (y_k1 - y_k3).__mul_neg_i();
        for(int k = 1; k < m; k++) {
            k0 = (k      )*s_out;
            k1 = (k +   m)*s_out;
            k2 = (k + 2*m)*s_out;
            k3 = (k + 3*m)*s_out;
            y_k0  = y[k0];
            y_k1  = y[k1];
            y_k2  = y[k2];
            y_k3  = y[k3];
            y[k0] = wk2.fmadd( y_k2, y_k0) + wk1.fmadd(y_k1, wk3*y_k3);
            y[k1] = wk2.fmadd(-y_k2, y_k0) + wk1.fmsub(y_k1, wk3*y_k3).__mul_i();
            y[k2] = wk2.fmadd( y_k2, y_k0) - wk1.fmadd(y_k1, wk3*y_k3);
            y[k3] = wk2.fmadd(-y_k2, y_k0) + wk3.fmsub(y_k3, wk1*y_k1).__mul_i();
            wk1 *= w1; wk2 *= w2; wk3 *= w3;
        }
    }
}

// External function to call the C-T radix-4 FFT
template<typename F, int L>
inline void pow4_FFT(Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, biFuncNode<F,L>* sRoot, direction dir) {
    const size_t N = sRoot->sz; // Size of problem
    pow4_FFT_helper(N, x, y, s_in, s_out, dir); // Call the radix-2 FFT
}

// External & Internal function for radix-N1 C-T FFTs
template<typename F, int L>
inline void composite_FFT(Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, biFuncNode<F,L>* sRoot, direction dir) {
    // Retrieve N
    size_t N = sRoot->sz;

    // Find the children on the call-graph
    biFuncNode<F,L>* left = sRoot + sRoot->left;
    biFuncNode<F,L>* right = sRoot + sRoot->right;

    // Get the size of the sub-problems
    size_t N1 = left->sz;
    size_t N2 = right->sz;

    // I'm currently using a temporary storage space malloc'd in recursive calls.
    // This isn't optimal and will change as the engine develops
    Complex<F,L>* z  = sRoot->workspace.data();
    // Find the FFT of the "rows" of the input signal and twiddle them accordingly
    Complex<F,L> w1 = omega<F,L>::get(1, N, dir);
    Complex<F,L> wj1 = Complex<F,L>(1., 0.);
    for(size_t j1 = 0; j1 < N1; j1++) {
        Complex<F,L> wk2 = wj1;
        right->fptr(&x[j1*s_in], &z[N2*j1], N1*s_in, 1, right, dir);
        if(j1 > 0) {
            for(size_t k2 = 1; k2 < N2; k2++) {
                z[j1*N2 + k2] *= wk2;
                wk2 *= wj1;
            }
        }
        wj1 *= w1;
    }

    /* Take the FFT of the "columns" of the output from the above "row" FFTs.
     * Don't need j1 for the second transform as it's just the indexer.
     * Take strides of N2 since z is allocated on the fly in this function for N.
     */
    for(size_t k2 = 0; k2 < N2; k2++) {
        left->fptr(&z[k2], &y[k2*s_out], N2, N2*s_out, left, dir);
    }
}

// A factoring function for the reference composite FFT
inline size_t referenceFactor(const size_t f) {
    // Check if it's even
    if((f & 0x1) == 0) return 2;

    // Check all odd numbers after that
    for(size_t k = 3; k*k < f; k+=2) {
        if( f % k == 0 ) return k;
    }

    // return f if no factor was found
    return f;
}

// Implementation for Rader's Algorithm
template<typename F, int L>
inline void rader_FFT(Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, biFuncNode<F,L>* sRoot, direction dir, size_t a, size_t ainv) {
    // Size of the problem
    size_t p = sRoot->sz;

    // Find the children on the call-graph
    biFuncNode<F,L>* subFFT = sRoot + sRoot->left;

    // Temporary workspace
    Complex<F,L>* z = sRoot->workspace.data();
    // Loop variables
    int ak = 1;
    int akinv = 1;
    Complex<F,L> y0 = x[0];

    // First, "invert" the order of x
    for(size_t k = 0; k < (p-1); k++) {
        y[k*s_out] = x[akinv*s_in];
        y0 = y0 + y[k*s_out];
        ak = (ak*a) % p;
        akinv = (akinv*ainv) % p;
    }

    // Convolve the resulting vector with twiddle vector

    // First fft the resulting shuffled vector
    subFFT->fptr(y, &z[0], s_out, 1, subFFT, dir);

    // Perform cyclic convolution
    for(size_t m = 0; m < (p-1); m++) {
        Complex<F,L> Cm = omega<F,L>::get(1, p, dir);
        ak = a;
        for(size_t k = 1; k < (p-1); k++) {
            Cm = Cm + omega<F,L>::get(p*(k*m+ak) - ak, p*(p-1), dir);
            ak = (ak*a) % p;
        }
        y[m*s_out] = z[m]*Cm;
    }

    // Bring back into signal domain
    subFFT->fptr(y, &z[0], s_out, 1, subFFT, (direction) (-1*((int) dir)));

    // Shuffle as needed
    ak = 1;
    y[0] = y0;
    for(size_t m = 0; m < (p-1); m++) {
        y[ak*s_out] = x[0] + (z[m]/((double) (p-1)));
        ak = (ak*a) % p;
    }
}

// Internal recursive helper-function that calculates the FFT of a signal with length 3^k
template<typename F, int L>
inline void pow3_FFT_helper(size_t N, Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, direction dir, Complex<F,L>& plus120, Complex<F,L>& minus120) {

    // Calculate the DFT manually if necessary
    if(N == 3) {
        y[0] = x[0] + x[s_in] + x[2*s_in];
        y[s_out] = x[0] + plus120*x[s_in] + minus120*x[2*s_in];
        y[2*s_out] = x[0] + minus120*x[s_in] + plus120*x[2*s_in];
        return;
    }

    // Calculate the size of the sub-problem
    size_t Nprime = N/3;

    // Divide into sub-problems
    pow3_FFT_helper(Nprime, x, y, s_in*3, s_out, dir, plus120, minus120);
    pow3_FFT_helper(Nprime, x+s_in, y+Nprime*s_out, s_in*3, s_out, dir, plus120, minus120);
    pow3_FFT_helper(Nprime, x+2*s_in, y+2*Nprime*s_out, s_in*3, s_out, dir, plus120, minus120);

    // Combine the sub-problem solutions
    Complex<F,L> w1 = omega<F,L>::get(1, N, dir);
    Complex<F,L> w2 = w1*w1;
    Complex<F,L> wk1 = w1;
    Complex<F,L> wk2 = w2;

    int k1 =                0;
    int k2 =   Nprime * s_out;
    int k3 = 2*Nprime * s_out;

    Complex<F,L> tmpk     = y[k1];
    Complex<F,L> tmpk_p_1 = y[k2];
    Complex<F,L> tmpk_p_2 = y[k3];

    y[k1] =                tmpk_p_2 +               tmpk_p_1+ tmpk;
    y[k2] = minus120.fmadd(tmpk_p_2, plus120.fmadd( tmpk_p_1, tmpk));
    y[k3] = plus120.fmadd( tmpk_p_2, minus120.fmadd(tmpk_p_1, tmpk));

    for(size_t k = 1; k < Nprime; k++) {
        // Index calculation
        k1 = k * s_out;
        k2 = (  Nprime + k) * s_out;
        k3 = (2*Nprime + k) * s_out;

        // Storing temporary variables
        tmpk     = y[k1];
        tmpk_p_1 = y[k2];
        tmpk_p_2 = y[k3];

        // Reassigning the output
        y[k1] = wk2.fmadd(           tmpk_p_2, wk1.fmadd(           tmpk_p_1, tmpk));
        y[k2] = wk2.fmadd(minus120 * tmpk_p_2, wk1.fmadd(plus120  * tmpk_p_1, tmpk));
        y[k3] = wk2.fmadd(plus120  * tmpk_p_2, wk1.fmadd(minus120 * tmpk_p_1, tmpk));

        // Twiddle factors
        wk1 *= w1; wk2 *= w2;
    }
}

// External-facing function for performing an FFT on signal with length N = 3^k
template<typename F, int L>
inline void pow3_FFT(Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, biFuncNode<F,L>* sRoot, direction dir) {
    const size_t N = sRoot->sz;
    Complex<F,L> plus120 (-0.5, -sqrt(3)/2.);
    Complex<F,L> minus120 (-0.5, sqrt(3)/2.);
    switch(dir) {
        case direction::forward:  pow3_FFT_helper(N, x, y, s_in, s_out, dir, plus120, minus120); break;
        case direction::backward: pow3_FFT_helper(N, x, y, s_in, s_out, dir, minus120, plus120); break;
    }
}

}
}

#endif // END HEFFTE_STOCK_ALGOS_H
