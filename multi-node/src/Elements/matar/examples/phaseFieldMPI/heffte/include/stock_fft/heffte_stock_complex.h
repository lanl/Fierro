/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_STOCK_COMPLEX_H
#define HEFFTE_STOCK_COMPLEX_H

#include <type_traits>
#include <iostream>

#include "heffte_stock_vec_types.h"

/**
 * This file holds the tools of using Complex numbers in a templated fashion
 * using vectorized intrinsics.
 */
namespace heffte {
namespace stock {
/*! \ingroup stockbackend
 * \brief Custom complex type taking advantage of vectorization
 * A Complex Type intrinsic to HeFFTe that takes advantage of vectorization.
 * It is currently compatible with AVX2 to create 1 double-precision complex,
 * 2 double-precision complex, 2 single-precision complex, or 4 single-precision
 * complex numbers.
 */
template<typename F, int L>
class alignas(L*sizeof(F)) Complex {
    public:
        // One 64-bit Complex-- 2 doubles-- pack<double, 2>::type == _m128d
        // Two 64-bit Complex-- 4 doubles-- pack<double, 4>::type == _m256d
        // Two 64-bit Complex-- 4 floats -- pack<float, 4>::type == _m128
        // Four 64-bit Complex-- 8 floats -- pack<float, 8>::type == _m256

        // Constructors for the Complex class
        //! \brief Load from an array of primitives
        explicit Complex(F* const f): var(mm_load<F,L>(f)) {}
        //! \brief Load from an initializer list of primitives
        explicit Complex(std::initializer_list<F> il): var(mm_load<F,L>(il.begin())) {};
        //! \brief Load from an existing vector pack
        explicit Complex(typename pack<F,L>::type v): var(v) {}
        //! \brief Load from real and imaginary parts (repeating for all numbers in pack)
        explicit Complex(F x, F y): var(mm_pair_set<F,L>(x, y)) {}
        //! \brief Load from pointer to existing std::complex numbers
        explicit Complex(std::complex<F> const *c): var(mm_complex_load<F,L>(c)) {}
        //! \brief Load from strided pointer to existing std::complex numbers
        explicit Complex(std::complex<F> const *c, int stride): var(mm_complex_load<F,L>(c, stride)) {}
        //! \brief Load from initializer list of existing std::complex numbers
        explicit Complex(std::initializer_list<std::complex<F>> il): var(il.size() == 1 ? mm_pair_set<F,L>((*il.begin()).real(), (*il.begin()).imag()) : mm_complex_load<F,L>(il.begin())) {};
        //! \brief Default constructor of zeros
        explicit Complex(): var(mm_zero<F,L>()) {}

        ///////////////////////////////////////////////////////////
        /* Basic operations with another pack of complex numbers */
        ///////////////////////////////////////////////////////////

        //! \brief Add with another pack of complex number
        Complex<F,L> operator+(Complex<F,L> const &o) {
            return Complex(mm_add(var, o.var));
        }

        //! \brief Subtract another pack of complex number
        Complex<F,L> operator-(Complex<F,L> const &o) {
            return Complex(mm_sub(var, o.var));
        }

        //! \brief Multiply by another pack of complex number
        Complex<F,L> operator*(Complex<F,L> const & o) {
            return Complex(mm_complex_mul(var, o.var));
        }

        //! \brief Divide by another pack of complex number
        Complex<F,L> operator/(Complex<F,L> const & o) {
            return Complex(mm_complex_div(var, o.var));
        }

        //! \brief Add with and assign another complex number
        Complex<F,L> operator+=(Complex<F,L> const &o) {
            var = mm_add(var, o.var);
            return *this;
        }

        //! \brief Subtract and assign another complex number from this
        Complex<F,L> operator-=(Complex<F,L> const &o) {
            var = mm_sub(var, o.var);
            return *this;
        }

        //! \brief Multiply by and assign another complex number
        Complex<F,L> operator*=(Complex<F,L> const &o) {
            var = mm_complex_mul(var, o.var);
            return *this;
        }

        //! \brief Divide by and assign another complex number
        Complex<F,L> operator/=(Complex<F,L> const &o) {
            var = mm_complex_div(var, o.var);
            return *this;
        }

        ///////////////////////////////////////////////////////////////
        /* Basic operations with a single real floating point number */
        ///////////////////////////////////////////////////////////////

        //! \brief Add with a floating point number
        Complex<F,L> operator+(F o) {
            return Complex(mm_add(var, mm_pair_set<F,L>(o, 0)));
        }

        //! \brief Subtract a floating point number
        Complex<F,L> operator-(F o) {
            return Complex(mm_sub(var, mm_pair_set<F,L>(o, 0)));
        }

        //! \brief Multiply by a floating point number
        Complex<F,L> operator*(F o) {
            return Complex(mm_mul(var, mm_set1<F,L>(o)));
        }

        //! \brief Divide by a floating point number
        Complex<F,L> operator/(F o) {
            return Complex(mm_div(var, mm_set1<F,L>(o)));
        }

        //! \brief Add with and assign a floating point number
        Complex<F,L> operator+=(F o) {
            var = mm_add(var, mm_pair_set<F,L>(o, 0));
            return *this;
        }

        //! \brief Subtract and assign a floating point number
        Complex<F,L> operator-=(F o) {
            var = mm_sub(var, mm_pair_set<F,L>(o, 0));
            return *this;
        }

        //! \brief Multiply by and assign a floating point number
        Complex<F,L> operator*=(F o) {
            var = mm_mul(var, mm_set1<F,L>(o));
            return *this;
        }

        // Divide by and assign a floating point number
        Complex<F,L> operator/=(F o) {
            var = mm_div(var, mm_set1<F,L>(o));
            return *this;
        }

        ///////////////////
        /* Other methods */
        ///////////////////

        //! \brief Fused multiply add
        Complex<F,L> fmadd(Complex<F,L> const & y, Complex<F,L> const & z) {
            return Complex(mm_complex_fmadd(var, y.var, z.var));
        }

        //! \brief Fused multiply add
        Complex<F,L> fmsub(Complex<F,L> const & y, Complex<F,L> const & z) {
            return Complex(mm_complex_fmsub(var, y.var, z.var));
        }

        //! \brief Negate the complex number
        Complex<F,L> operator-() {
            return Complex<F,L>(mm_neg(var));
        }

        //! \brief Store the modulus of the complex number in an array of size L/2
        void modulus(F* dest) {
            typename pack<F, L>::type res = mm_complex_mod(var);
            for(int i = 0; i < L/2; i++) {
                dest[i] = res[i*2];
            }
        }

        //! \brief Return the modulus of the complex number in a vector pack
        typename pack<F,L>::type modulus() {
            return mm_complex_mod(var);
        }

        //! \brief Conjugate the current complex number
        Complex<F,L> conjugate() {
            return Complex(mm_complex_conj(var));
        }

        //! \brief Multiply the complex number by i
        Complex<F,L> __mul_i() {
            return Complex(mm_complex_mul_i(var));
        }

        //! \brief Multiply the complex number by i
        Complex<F,L> __mul_neg_i() {
            return Complex(mm_complex_mul_neg_i(var));
        }

        //! \brief Store the Complex number in an array of length L
        void get(F* dest) {
            mm_store<F,L>(dest, var);
        }

        //! \brief Access the ith Complex number as a std::complex
        std::complex<F> operator[](std::size_t idx) {
            return std::complex<F>(reinterpret_cast<F*>(&var)[2*idx], reinterpret_cast<F*>(&var)[2*idx + 1]);
        }

        //! \brief Return a vector pack representation of this number
        typename pack<F,L>::type get() const {
            return var;
        }

    private:
        //! \brief Vector pack representing the complex number(s)
        typename pack<F,L>::type var;
};




//! \brief Determines how the complex number should be printed to an ostream
template<typename F, int L>
inline std::ostream& operator<<(std::ostream& os, const Complex<F,L>& dt){
    typename pack<F, L>::type var = dt.get();
    os << "( ";
    for(int i = 0; i < L; i+=2) {
        os << var[i];
        if(var[i+1] < 0) os << " - " << -var[i+1] << "i";
        else os << " + " << var[i+1] << "i";
        if(i+2 < L) os << ", ";
    }
    os << " )";
    return os;
}

}
}

#endif // HEFFTE_STOCK_COMPLEX_H
