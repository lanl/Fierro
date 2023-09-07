/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_STOCK_TREE_H
#define HEFFTE_STOCK_TREE_H

#include "heffte_stock_algos.h"
#include <iostream>

//! \brief This is the minimum signal size for us to use Rader's algorithm
#define HEFFTE_STOCK_RADER_MIN 100000

namespace heffte {
namespace stock {

//! \brief Const array of (NOT NECESSARILY PRIME) factors
static const std::array<size_t,200> factors {   4,    2,    3,    5,    7,    11,   13,   16,   17,   19,
                                                23,   29,   31,   37,   41,   43,   47,   53,   59,   61,
                                                67,   71,   73,   79,   83,   89,   97,   101,  103,  107,
                                                109,  113,  127,  131,  137,  139,  149,  151,  157,  163,
                                                167,  173,  179,  181,  191,  193,  197,  199,  211,  223,
                                                227,  229,  233,  239,  241,  251,  257,  263,  269,  271,
                                                277,  281,  283,  293,  307,  311,  313,  317,  331,  337,
                                                347,  349,  353,  359,  367,  373,  379,  383,  389,  397,
                                                401,  409,  419,  421,  431,  433,  439,  443,  449,  457,
                                                461,  463,  467,  479,  487,  491,  499,  503,  509,  521,
                                                523,  541,  547,  557,  563,  569,  571,  577,  587,  593,
                                                599,  601,  607,  613,  617,  619,  631,  641,  643,  647,
                                                653,  659,  661,  673,  677,  683,  691,  701,  709,  719,
                                                727,  733,  739,  743,  751,  757,  761,  769,  773,  787,
                                                797,  809,  811,  821,  823,  827,  829,  839,  853,  857,
                                                859,  863,  877,  881,  883,  887,  907,  911,  919,  929,
                                                937,  941,  947,  953,  967,  971,  977,  983,  991,  997,
                                                1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061,
                                                1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123,
                                                1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213};

//! \brief Find the smallest usable factor of f (could be done at compile-time)
inline size_t factor(const size_t f) {
    if(f < HEFFTE_STOCK_THRESHOLD) return f;
    size_t k = 0;
    // Prioritize factors in the factors array
    for(; k < factors.size(); k++) {
        if( f % factors[k] == 0) return factors[k];
    }

    // If none of those work, turn to odd numbers that are greater than the last index
    for(k = factors[k - 1]; k*k < f; k+=2) {
        if( f % k == 0 ) return k;
    }

    // return f if no factor was found
    return f;
}

//! \brief Find all the factors to test for primitive root
inline std::vector<size_t> findFactorRader(size_t factor) {
    size_t f = factor;
    std::vector<size_t> ret {};
    if((f % 2) == 0) {
        while((f % 2) == 0) {
            f /= 2;
        }
        ret.push_back(factor/2);
    }

    for(size_t k = 3; k <= f; k += 2) {
        if((f % k) == 0) {
            while((f % k) == 0) f /= k;
            ret.push_back(factor/k);
        }
    }
    return ret;
}

//! \brief Get (base ^ pow) % mod (could be done at compile-time)
inline size_t modPow(size_t base, size_t pow, size_t mod) {
    size_t ret = 1;
    while(pow > 0) {
        if((pow & 0x1) != 0) {
            ret = (ret*base) % mod;
        }
        base = (base*base) % mod;
        pow >>= 1;
    }
    return ret;
}

//! \brief Find the prime root for p, prime
inline size_t primeRoot(size_t p) {
    size_t phi = p - 1;
    std::vector<size_t> f_vec = findFactorRader(phi);
    bool ret = false;
    size_t curr = 1;
    while(!ret) {
        curr++;
        for(auto& i : f_vec) {
            size_t pow = modPow(curr, i, p);
            if(pow == 1) {
                ret = false;
                break;
            }
            ret = true;
        }
    }
    return curr;
}

//! \brief Check if n is a power of pow (could be done at compile-time)
inline bool power_of(const size_t n, const size_t pow) {
    if(n == 1) return false;
    size_t k = n;
    if(pow == 2) {
        unsigned char sum = 0x1 & n;
        for(;k;k>>=1) sum += 0x1&k;
        return sum == 1;
    }
    while(!(k % pow)) k/= pow;
    return k == 1;
}

/*! \ingroup stockbackend
 *  \brief Count nodes in a call-graph for FFT of length N (could be done at compile-time)
 *  First check if N is a power of something. Then, check if it is divisible
 *  by a power, and finally simply just factor.
 */
inline size_t numNodesFactorHelper(const size_t N) {
    // Check if power
    if(power_of(N, 4) ||
        power_of(N, 2) ||
        power_of(N, 3)) {
        return N;
    }
    // Check if divisible by a power
    if((N % 4) == 0) {
        size_t k = N;
        while ((k % 4) == 0) k /= 4;
        return N / k;
    }
    if((N % 2) == 0) {
        size_t k = N;
        while ((k % 2) == 0) k /= 2;
        return N / k;
    }
    if((N % 3) == 0) {
        size_t k = N;
        while ((k % 3) == 0) k /= 3;
        return N / k;
    }
    // Take a factor
    size_t k = factor(N);
    return (k == N && k > HEFFTE_STOCK_RADER_MIN) ? N-1 : k;
}

//! \brief Divide N by k if N composite, otherwise return zero (could be done at compile-time)
inline size_t getLeftover(const size_t N, const size_t k) {
    return (k == N-1) ? 0 : N/k;
}

//! \brief Return type of FFT and initialize factor in composite case. See NumNodesHelper for info.
inline std::pair<fft_type, size_t> fptrFactorHelper(const size_t N) {
    if(N < HEFFTE_STOCK_THRESHOLD) {
        return std::pair<fft_type,size_t> {fft_type::discrete, N};
    }
    // Check if N is a power
    if(power_of(N, 4)) {
        return std::pair<fft_type,size_t> {fft_type::pow4, N};
    }
    if(power_of(N, 2)) {
        return std::pair<fft_type,size_t> {fft_type::pow2, N};
    }
    if(power_of(N, 3)) {
        return std::pair<fft_type,size_t> {fft_type::pow3, N};
    }
    // Check if N is divisible by a power
    if((N % 4) == 0) {
        size_t ell = N;
        while ((ell % 4) == 0) ell /= 4;
        return std::pair<fft_type,size_t> {fft_type::composite, N/ell};
    }
    if((N % 2) == 0) {
        size_t ell = N;
        while ((ell % 2) == 0) ell /= 2;
        return std::pair<fft_type,size_t> {fft_type::composite, N/ell};
    }
    if((N % 3) == 0) {
        size_t ell = N;
        while ((ell % 3) == 0) ell /= 3;
        return std::pair<fft_type,size_t> {fft_type::composite, N/ell};
    }
    // Factor N
    size_t k = factor(N);
    if(k == N) {
        // Worry about whether to call Rader's algorithm
        if(N > HEFFTE_STOCK_RADER_MIN) {
            k = N-1;
            return std::pair<fft_type,size_t> {fft_type::rader, k};
        }
        return std::pair<fft_type,size_t> {fft_type::discrete, k};
    }
    return std::pair<fft_type,size_t> {fft_type::composite, k};
}

//! \brief Recursively initialize a call-graph given an array (could be done at compile-time)
template<typename F, int L>
inline size_t init_fft_tree(biFuncNode<F,L>* sRoot, const size_t N) {
    std::pair<fft_type, size_t> pr = fptrFactorHelper(N);
    fft_type type = pr.first; size_t k = pr.second;
    if(type == fft_type::rader) {
        size_t a = primeRoot(N);
        size_t ainv = modPow(a, N-2, N);
        *sRoot = biFuncNode<F,L>(N, a, ainv);
    }
    else {
        *sRoot = biFuncNode<F,L>(N, type);
    }
    if(type == fft_type::discrete ||
        type == fft_type::pow2     ||
        type == fft_type::pow3     ||
        type == fft_type::pow4) {
        return 1;
    }
    size_t q = getLeftover(N, k);
    size_t l = init_fft_tree(sRoot + 1, k);
    size_t r = (type == fft_type::rader) ? 0 : init_fft_tree(sRoot + 1 + l, q);
    sRoot->left = 1;
    sRoot->right = 1 + l;
    return 1 + l + r;
}

//! \brief Get the number of nodes appropriate for a given length N (could be done at compile-time)
inline size_t getNumNodes(const size_t N) {
    size_t k = numNodesFactorHelper(N);
    size_t rem = getLeftover(N, k);
    if(k == N) return 1;
    size_t left_nodes = getNumNodes(k);
    size_t right_nodes = (rem == 0) ? 0 : getNumNodes(rem);
    return 1 + left_nodes + right_nodes;
}

}
}

#endif // END HEFFTE_STOCK_TREE_H