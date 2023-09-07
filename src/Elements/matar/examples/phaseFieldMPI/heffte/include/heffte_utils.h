/*
    -- heFFTe --
       Univ. of Tennessee, Knoxville
       @date
*/

#ifndef HEFFTE_UTILS_H
#define HEFFTE_UTILS_H

#include <algorithm>
#include <vector>
#include <complex>
#include <memory>
#include <numeric>
#include <algorithm>
#include <functional>
#include <cassert>
#include <utility>
#include <iostream>
#include <ostream>
#include <iomanip>
#include <string>
#include <deque>
#include <fstream>
#include <mpi.h>

#include "heffte_config.h"

namespace heffte {

/*!
 * \ingroup fft3dmisc
 * \brief Replace with the C++ 2014 std::exchange later.
 */
template<class T, class U = T>
T c11_exchange(T& obj, U&& new_value)
{
    T old_value = std::move(obj);
    obj = std::forward<U>(new_value);
    return old_value;
}

/*!
 * \ingroup fft3d
 * \addtogroup hefftempi Helper wrappers around MPI methods
 *
 * HeFFTe is using the C-style of API for the message passing interface (MPI).
 * A set of inline wrappers allow for easier inclusion into C++ methods
 * which helps preserve const-correctness and simplifies variable initialization.
 */

/*!
 * \ingroup hefftempi
 * \brief Wrappers to miscellaneous MPI methods giving a more C++-ish interface.
 */
namespace mpi {

/*!
 * \ingroup hefftempi
 * \brief Returns the rank of this process within the specified \b comm.
 *
 * \param comm is an MPI communicator associated with the process.
 *
 * \returns the rank of the process within the \b comm
 *
 * Uses MPI_Comm_rank().
 *
 * Example:
 * \code
 *  int const me = mpi::comm_rank(comm); // good C++ style
 * \endcode
 * as opposed to:
 * \code
 *  int me; // uninitialized variable, cannot be const
 *  MPI_Comm_rank(comm, &me); // initialization takes a second line, bad C++ style
 * \endcode
 */
inline int comm_rank(MPI_Comm const comm){
    int me;
    MPI_Comm_rank(comm, &me);
    return me;
}
/*!
 * \ingroup hefftempi
 * \brief Returns \b true if this process has the \b me rank within the MPI_COMM_WORLD (useful for debugging).
 *
 * \param rank within the world communicator
 * \return true if the rank of this processor matches the given \b rank
 *
 * Useful for debugging, e.g., cout statements can be easily constrained to a single processor, e.g.,
 * \code
 *  if (mpi::world_rank(3)){
 *      cout << ....; // something about the state
 *  }
 * \endcode
 */
inline bool world_rank(int rank){ return (comm_rank(MPI_COMM_WORLD) == rank); }
/*!
 * \ingroup hefftempi
 * \brief Returns the rank of this process within the MPI_COMM_WORLD (useful for debugging).
 *
 * Useful for debugging.
 */
inline int world_rank(){ return comm_rank(MPI_COMM_WORLD); }

/*!
 * \ingroup hefftempi
 * \brief Write the message and the data from the vector-like \b x, performed only on rank \b me (if positive), otherwise using all ranks.
 *
 * Very useful for debugging purposes, when a vector or vector-like object has to be inspected for a single MPI rank.
 * \tparam vector_like is an object that can be used for ranged for-loop, e.g., std::vector or std::array
 * \param me the rank to write to cout
 * \param x is the data to be written out
 * \param message will be written on the line before \b x, helps identify what \b x should contain, e.g., result or reference data
 */
template<typename vector_like>
void dump(int me, vector_like const &x, std::string const &message){
    if (me < 0 or world_rank(me)){
        std::cout << message << "\n";
        for(auto i : x) std::cout << i << "  ";
        std::cout << std::endl;
    }
}

/*!
 * \ingroup hefftempi
 * \brief Returns the size of the specified communicator.
 *
 * \param comm is an MPI communicator associated with the process.
 *
 * \returns the number of ranks associated with the communicator (i.e., size).
 *
 * Uses MPI_Comm_size().
 */
inline int comm_size(MPI_Comm const comm){
    int nprocs;
    MPI_Comm_size(comm, &nprocs);
    return nprocs;
}
/*!
 * \ingroup hefftempi
 * \brief Creates a new sub-communicator from the provided processes in \b comm.
 *
 * \param ranks is a list of ranks associated with the \b comm.
 * \param comm is an active communicator that holds all processes in the \b ranks.
 *
 * \returns a new communicator that uses the selected ranks.
 *
 * Uses MPI_Comm_group(), MPI_Group_incl(), MPI_Comm_create(), MPI_Group_free().
 */
inline MPI_Comm new_comm_from_group(std::vector<int> const &ranks, MPI_Comm const comm){
    MPI_Group orig_group, new_group;
    MPI_Comm_group(comm, &orig_group);
    MPI_Group_incl(orig_group, (int) ranks.size(), ranks.data(), &new_group);
    MPI_Comm result;
    MPI_Comm_create(comm, new_group, &result);
    MPI_Group_free(&orig_group);
    MPI_Group_free(&new_group);
    return result;
}

/*!
 * \ingroup hefftempi
 * \brief Calls free on the MPI comm.
 *
 * \param comm is the communicator to be deleted, cannot be used after this call.
 *
 * Uses MPI_Comm_free().
 *
 * Note that the method would use const_cast() to pass the MPI_Comm (which is a pointer)
 * to the delete method. This circumvents the C-style of API that doesn't respect the fact
 * that deleting a const-pointer is an acceptable operation.
 */
inline void comm_free(MPI_Comm const comm){
    if (MPI_Comm_free(const_cast<MPI_Comm*>(&comm)) != MPI_SUCCESS)
        throw std::runtime_error("Could not free a communicator.");
}

/*!
 * \ingroup hefftempi
 * \brief Returns the MPI equivalent of the \b scalar C++ type.
 *
 * This template cannot be instantiated directly, only the specializations are useful.
 * Direct instantiation indicates an unknown conversion from a C++ to an MPI type.
 * \tparam scalar a C++ scalar type, e.g., float, double, std::complex<float>, etc.
 *
 * \returns the MPI equivalent, e.g., MPI_FLOAT, MPI_DOUBLE, MPI_C_COMPLEX, etc.
 */
template<typename scalar> inline MPI_Datatype type_from(){
    // note that "!std::is_same<scalar, scalar>::value" is always false,
    // but will not be checked until the template is instantiated
    static_assert(!std::is_same<scalar, scalar>::value, "The C++ type has unknown MPI equivalent.");
    return MPI_BYTE; // some compilers complain about lack of return statement.
}
/*!
 * \ingroup hefftempi
 * \brief Specialization to hand the int type.
 */
template<> inline MPI_Datatype type_from<int>(){ return MPI_INT; }
/*!
 * \ingroup hefftempi
 * \brief Specialization to hand the float type.
 */
template<> inline MPI_Datatype type_from<float>(){ return MPI_FLOAT; }
/*!
 * \ingroup hefftempi
 * \brief Specialization to hand the double type.
 */
template<> inline MPI_Datatype type_from<double>(){ return MPI_DOUBLE; }
/*!
 * \ingroup hefftempi
 * \brief Specialization to hand the single-precision complex type.
 */
template<> inline MPI_Datatype type_from<std::complex<float>>(){ return MPI_C_COMPLEX; }
/*!
 * \ingroup hefftempi
 * \brief Specialization to hand the double-precision complex type.
 */
template<> inline MPI_Datatype type_from<std::complex<double>>(){ return MPI_C_DOUBLE_COMPLEX; }

}

/*!
 * \ingroup fft3d
 * \brief Struct to specialize to allow HeFFTe to recognize custom single precision complex types.
 *
 * Specializations of this struct will allow HeFFTe to recognize custom complex types
 * that are ABI compatible with std::complex.
 * In this context, ABI compatible means that it is safe to use reinterpret_cast
 * between raw-arrays arrays of the two types.
 *
 * \tparam scalar_type indicates the type in question, if the type is ABI compatible
 *          with single precision complex std::complex<float> then the specialization
 *          must inherit from std::true_type, otherwise inherit from std::false_type.
 *          Note that the true/false types define a static const bool member
 *          called value that is correspondingly true/false.
 *
 * See std::is_zcomplex for specialization for double precision complex types,
 * the ccomplex and zcomplex names are mimicking by the BLAS naming conventions, e.g., cgemm() and zgemm().
 *
 * Example:
 * \code
 *  struct custom_single_precision_complex{
 *      float real, imag;
 *  }
 *  namespace heffte {
 *      template<> struct is_ccomplex<custom_single_precision_complex> : std::true_type{};
 *  }
 * \endcode
 */
template<typename scalar_type> struct is_ccomplex : std::false_type{};
/*!
 * \ingroup fft3d
 * \brief Struct to specialize to allow HeFFTe to recognize custom double precision complex types.
 *
 * Specializations of this struct will allow HeFFTe to recognize custom complex types
 * that are ABI compatible with std::complex<double>. See heffte::is_ccomplex for details.
 *
 * Example:
 * \code
 *  struct custom_double_precision_complex{
 *      double real, imag;
 *  }
 *  namespace heffte {
 *      template<> struct is_ccomplex<custom_double_precision_complex> : std::true_type{};
 *  }
 * \endcode
 */
template<typename scalar_type> struct is_zcomplex : std::false_type{};

/*!
 * \ingroup fft3dcomplex
 * \brief By default, HeFFTe recognizes std::complex<float>.
 */
template<> struct is_ccomplex<std::complex<float>> : std::true_type{};
/*!
 * \ingroup fft3dcomplex
 * \brief By default, HeFFTe recognizes std::complex<double>.
 */
template<> struct is_zcomplex<std::complex<double>> : std::true_type{};

/*!
 * \ingroup fft3dcomplex
 * \brief Struct to specialize that returns the C++ equivalent of each type.
 *
 * Given a type that is either float, double, or recognized by the heffte::is_ccomplex
 * and heffte::is_zcomplex templates, this struct will define a member type (called type)
 * that will define the corresponding C++ equivalent.
 *
 * Example:
 * \code
 *  struct custom_double_precision_complex{
 *      double real, imag;
 *  }
 *  namespace heffte {
 *      template<> struct is_ccomplex<custom_double_precision_complex> : std::true_type{};
 *  }
 *  ...
 *  static_assert(std::is_same<typename define_standard_type<custom_double_precision_complex>::type,
 *                             std::complex<double>>::value,
 *                "error: custom_double_precision_complex not equivalent to std::complex<double>");
 *
 *  template<typename input_type>
 *  void foo(input_type x[]){
 *      auto y = reinterpret_cast<typename define_standard_type<input_type>::type*>(x);
 *      ...
 *  }
 *  ...
 *  std::vector<custom_double_precision_complex> x(10);
 *  foo(x.data()); // here input_type will be deduced to custom_double_precision_complex
 *                 // and inside foo() y will be std::complex<double>*
 * \endcode
 */
template<typename, typename = void> struct define_standard_type{};

/*!
 * \ingroup fft3dcomplex
 * \brief Type float is equivalent to float.
 */
template<> struct define_standard_type<float, void>{
    //! \brief Float is equivalent to float.
    using type = float;
};
/*!
 * \ingroup fft3dcomplex
 * \brief Type double is equivalent to double.
 */
template<> struct define_standard_type<double, void>{
    //! \brief Double is equivalent to double.
    using type = double;
};

/*!
 * \ingroup fft3dcomplex
 * \brief Every type with specialization of heffte::is_ccomplex to std::true_type is equivalent to std::complex<float>.
 */
template<typename scalar_type> struct define_standard_type<scalar_type, typename std::enable_if<is_ccomplex<scalar_type>::value>::type>{
    //! \brief If heffte::is_ccomplex is true_type, then the type is equivalent to std::complex<float>.
    using type =  std::complex<float>;
};
/*!
 * \ingroup fft3dcomplex
 * \brief Every type with specialization of heffte::is_zcomplex to std::true_type is equivalent to std::complex<double>.
 */
template<typename scalar_type> struct define_standard_type<scalar_type, typename std::enable_if<is_zcomplex<scalar_type>::value>::type>{
    //! \brief If heffte::is_ccomplex is true_type, then the type is equivalent to std::complex<double>.
    using type =  std::complex<double>;
};

/*!
 * \brief Converts an array of some type to an array of the C++ equivalent type.
 */
template<typename scalar_type>
typename define_standard_type<scalar_type>::type* convert_to_standard(scalar_type input[]){
    return reinterpret_cast<typename define_standard_type<scalar_type>::type*>(input);
}
/*!
 * \ingroup fft3dcomplex
 * \brief Converts a const array of some type to a const array of the C++ equivalent type.
 */
template<typename scalar_type>
typename define_standard_type<scalar_type>::type const* convert_to_standard(scalar_type const input[]){
    return reinterpret_cast<typename define_standard_type<scalar_type>::type const*>(input);
}

/*!
 * \ingroup fft3d
 * \addtogroup fft3dmisc Miscellaneous helpers
 *
 * Simple helper templates.
 */

/*!
 * \ingroup fft3dmisc
 * \brief Return the index of the last active (non-null) unique_ptr.
 *
 * The method returns -1 if all shapers are null.
 */
template<typename some_class>
int get_last_active(std::array<std::unique_ptr<some_class>, 4> const &shaper){
    int last = -1;
    for(int i=0; i<4; i++) if (shaper[i]) last = i;
    return last;
}

/*!
 * \ingroup fft3dmisc
 * \brief Return the number of active (non-null) unique_ptr.
 */
template<typename some_class>
int count_active(std::array<std::unique_ptr<some_class>, 4> const &shaper){
    int num = 0;
    for(int i=0; i<4; i++) if (shaper[i]) num++;
    return num;
}

/*!
 * \ingroup fft3dmisc
 * \brief Returns the max of the box_size() for each of the executors.
 */
template<typename some_class>
size_t get_max_box_size(std::array<some_class, 3> const &executors){
    size_t max_size = (executors[0]) ? executors[0]->box_size() : 0;
    max_size = std::max(max_size, (executors[1]) ? executors[1]->box_size() : static_cast<size_t>(0));
    return std::max(max_size, (executors[2]) ? executors[2]->box_size() : static_cast<size_t>(0));
}

/*!
 * \ingroup fft3dmisc
 * \brief Returns the max of the box_size() for each of the executors.
 */
template<typename some_class>
size_t get_max_box_size_r2c(std::array<some_class, 3> const &executors){
    size_t max_size = (executors[0]) ? executors[0]->complex_size() : 0;
    max_size = std::max(max_size, (executors[1]) ? executors[1]->box_size() : static_cast<size_t>(0));
    return std::max(max_size, (executors[2]) ? executors[2]->box_size() : static_cast<size_t>(0));
}
/*!
 * \ingroup fft3dmisc
 * \brief Returns the max of the workspace_size() for each of the executors.
 */
template<typename some_class>
size_t get_max_work_size(std::array<some_class, 3> const &executors){
    size_t max_size = (executors[0]) ? executors[0]->workspace_size() : 0;
    max_size = std::max(max_size, (executors[1]) ? executors[1]->workspace_size() : static_cast<size_t>(0));
    return std::max(max_size, (executors[2]) ? executors[2]->workspace_size() : static_cast<size_t>(0));
}

/*!
 * \ingroup fft3dmisc
 * \brief Returns the max of the workspace_size() for each of the executors.
 */
template<typename some_class_r2c, typename some_class>
size_t get_max_work_size(some_class_r2c const &executors_r2c, std::array<some_class, 2> const &executors){
    return std::max(executors_r2c->workspace_size(), std::max(executors[0]->workspace_size(), executors[1]->workspace_size()));
}

}

#endif /* HEFFTE_UTILS_H */
