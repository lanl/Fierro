
#include "swage.h"


//using namespace ViewCArrayMath;


// --------------------------------
// Math functors for Arrays
// --------------------------------

// forward declare classes

// CArray math
namespace ViewCArrayMath {
    
    // Functor Add(A, B), adding two 2D arrays
    class Add22 {
        
        ViewCArray <double> A_;
        ViewCArray <double> B_;
        
    public:
        
        // default constructor
        Add22();
        
        // constructor
        Add22(ViewCArray <double> &A, ViewCArray <double> &B);
        
        // calculate C=A+B via the interface A2D_add(C)
        void operator()(ViewCArray <double> &C_);
        
    }; // end functor
    
    
    // Functor Multiply11(a, b) where arrays a and b are 1D array
    class Multiply11 {
        
        ViewCArray <double> a_;
        ViewCArray <double> b_;
        
    public:
        
        // default constructor
        Multiply11();
        
        // constructor
        Multiply11(ViewCArray <double> &a, ViewCArray <double> &b);
        
        // calculate c=a*b via the interface A1D_dot_A1D(C)
        void operator()(double &c_);
        
    };
    
    // Functor Multiply21(A, b) where array A is 2D and b is a 1D array
    class Multiply21 {
        
        ViewCArray <double> A_;
        ViewCArray <double> b_;
        
    public:
        
        // default constructor
        Multiply21();
        
        // constructor
        Multiply21(ViewCArray <double> &A, ViewCArray <double> &b);
        
        // calculate c=A*b via the interface A2D_dot_A1D(c)
        void operator()(ViewCArray <double> &c_);
        
    };
    
    // Functor Multiply22(A, B) where arrays A and B are 2D
    class Multiply22 {
        
        ViewCArray <double> A_;
        ViewCArray <double> B_;
        
    public:
        
        // default constructor
        Multiply22();
        
        // constructor
        Multiply22(ViewCArray <double> &A, ViewCArray <double> &B);
        
        // calculate C=A*B via the interface A2D_dot_A2D(C)
        void operator()(ViewCArray <double> &C_);
        
    };
    
}; // end of namespace



// constructor overload
ViewCArrayMath::Add22::Add22(ViewCArray <double> &A, ViewCArray <double> &B){
    A_ = A;
    B_ = B;
}; // end constructor

void ViewCArrayMath::Add22::operator()(ViewCArray <double> &C_){
    
    assert(A_.order() == 2 && "ViewCArray Add22 dimensions do not match, first input array should be a 2D array!");
    assert(B_.order() == 2 && "ViewCArray Add22 dimensions do not match, second input array should be a 2D array!");
    assert(C_.order() == 2 && "ViewCArray Add22 dimensions do not match, second input array should be a 2D array!");
    assert(C_.dims(0) == A_.dims(0) && "ViewCArray Add22 dimensions do not match, output and first input array have different first dimensions!");
    assert(C_.dims(1) == B_.dims(1) && "ViewCArray Add22 dimensions do not match, output and second input array have different second dimensions!");
    assert(A_.dims(0) == B_.dims(0) && "ViewCArray Add22 dimensions do not match, input arrays A(i,j)+B(i,j) have wong i dimension!");
    assert(A_.dims(1) == B_.dims(1) && "ViewCArray Add22 dimensions do not match, input arrays A(i,j)+B(i,j) have wong j dimension!");
    
    // execute functor Multiply
    for (size_t i=0; i<A_.dims(0); i++){
        for (size_t j=0; j<A_.dims(1); j++){
            C_(i,j) = A_(i,j) + B_(i,j);
        }
    } // end for
}; // end functor overload


// constructor overload
ViewCArrayMath::Multiply11::Multiply11(ViewCArray <double> &a, ViewCArray <double> &b){
    a_ = a;
    b_ = b;
}; // end constructor

void ViewCArrayMath::Multiply11::operator()(double &c_){
    
    assert(a_.order() == 1 && "ViewCArray Multiply11 dimensions do not match, first input array should be a 1D array!");
    assert(b_.order() == 1 && "ViewCArray Multiply11 dimensions do not match, second input array should be a 1D array!");
    
    // initialize to zero
    c_ = 0.0;
    
    // execute functor Multiply
    for (size_t i=0; i<a_.dims(0); i++){
        c_ += a_(i) * b_(i);
    } // end for
}; // end functor overload

// constructor overload
ViewCArrayMath::Multiply21::Multiply21(ViewCArray <double> &A, ViewCArray <double> &b){
    A_ = A;
    b_ = b;
}; // end constructor


void ViewCArrayMath::Multiply21::operator()(ViewCArray <double> &c_){
    
    assert(A_.order() == 2 && "ViewCArray Multiply21 dimensions do not match, first input array should be a 2D array!");
    assert(b_.order() == 1 && "ViewCArray Multiply21 dimensions do not match, second input array should be a 1D array!");
    assert(c_.order() == 1 && "ViewCArray Multiply21 dimensions do not match, second input array should be a 1D array!");
    assert(c_.dims(0) == A_.dims(0) && "ViewCArray Multiply21 dimensions do not match, output and first input array have different first dimensions!");
    assert(c_.dims(1) == 0 && "ViewCArray Multiply21 dimensions do not match, output array should be a 1D array!");
    assert(A_.dims(1) == b_.dims(0) && "ViewCArray Multiply21 dimensions do not match, input arrays A(i,j)*b(j) have wong j dimension!");
    
    // initialize to zero
    for (size_t i=0; i<c_.dims(0); i++){
        c_(i) = 0.0;
    } // end for
    
    // execute functor Multiply
    for (size_t i=0; i<A_.dims(0); i++){
        for (size_t j=0; j<A_.dims(1); j++){
            c_(i) += A_(i,j) * b_(j);
        }
    } // end for
}; // end functor overload


// constructor overload
ViewCArrayMath::Multiply22::Multiply22(ViewCArray <double> &A, ViewCArray <double> &B){
    A_ = A;
    B_ = B;
}; // end constructor

void ViewCArrayMath::Multiply22::operator()(ViewCArray <double> &C_){
    
    assert(A_.order() == 2 && "ViewCArray Multiply22 dimensions do not match, first input array should be a 2D array!");
    assert(B_.order() == 2 && "ViewCArray Multiply22 dimensions do not match, second input array should be a 2D array!");
    assert(C_.order() == 2 && "ViewCArray Multiply22 dimensions do not match, second input array should be a 2D array!");
    assert(C_.dims(0) == A_.dims(0) && "ViewCArray Multiply22 dimensions do not match, output and first input array have different first dimensions!");
    assert(C_.dims(1) == B_.dims(1) && "ViewCArray Multiply22 dimensions do not match, output and second input array have different second dimensions!");
    assert(A_.dims(1) == B_.dims(0) && "ViewCArray Multiply22 dimensions do not match, input arrays A(i,j)*B(j,k) have wong j dimension!");
    
    // initialize to zero
    for (size_t i=0; i<C_.dims(0); i++){
        for (size_t j=0; j<C_.dims(1); j++){
            C_(i,j) = 0.0;
        }
    } // end for
    
    // execute functor Multiply
    for (size_t i=0; i<A_.dims(0); i++){
        for (size_t k=0; k<B_.dims(1); k++){
            for (size_t j=0; j<A_.dims(1); j++){
                C_(i,k) += A_(i,j) * B_(j,k);
            }
        }
    } // end for
}; // end functor overload



/** Test slam math implementation */
int main() {
  
    // a trad. array for testing
    double A[9] = {1,1,1,1,1,1,1,1,1};  // A initialized to 1
    double B[9] = {1,1,1,1,1,1,1,1,1};  // B initialized to 1
    double C[9]; // send results to C as in C=A+B
    
    
    // view as 2D arrays
    auto A_array = ViewCArray <double> (A, 3, 3);
    auto B_array = ViewCArray <double> (B, 3, 3);
    auto C_array = ViewCArray <double> (C, 3, 3);
    
    
    // Calculate C=A+B
    C_array = ViewCArrayMath::Add22(A_array, B_array);
    std::cout << "testing add \n";
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            std::cout << C_array(i,j) << "\n";
        }
    }
    std::cout << " " << std::endl;
    
    
    // Calculate C=A*B
    C_array = ViewCArrayMath::Multiply22(A_array, B_array);
    std::cout << "testing multiply \n";
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            std::cout << C_array(i,j) << "\n";
        }
    }
    
    std::cout << "--- finished ---" << std::endl;

  return 0;
}
