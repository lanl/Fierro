#include <stdio.h>
#include <array>
#include <variant>
#include <chrono>

#include "matar.h"

using namespace mtr; // matar namespace

// helper type for selecting variant type set by user
template<class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
// explicit deduction guide (not needed as of C++20)
template<class... Ts> overloaded(Ts...) -> overloaded<Ts...>;



// A notional class
class Data {
private:
    u_int nx_;
    u_int ny_;
    u_int nz_;
    
    CArrayKokkos <int> arr3D_;
    
public:
    
    // default constructor
    Data();
    
    // overload constructor to set dimensions
    Data(u_int nx, u_int ny, u_int nz);
    
    void some_fcn();
    
}; // end class Data

Data::Data(){};

Data::Data(u_int nx, u_int ny, u_int nz){
    
    nx_ = nx;
    ny_ = ny;
    nz_ = nz;
    
    arr3D_ = CArrayKokkos <int> (nx_, ny_, nz_);
};


void Data::some_fcn(){
    
    // parallel loop inside a class
    // The KOKKOS_CLASS_LAMBDA is [=, *this]. The *this in the lambda
    // capture gives access to the class data
    Kokkos::parallel_for("3DCArray",
                         Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0}, {nx_, ny_, nz_}),
                         KOKKOS_CLASS_LAMBDA(const int i, const int j, const int k) {
                             int idx = (i-1) * nx_ * ny_ + (j-1) * nz_ + (k-1);
                             arr3D_(i, j, k) = idx;
                         });
    Kokkos::fence();
    
    // now using the macros for a parallel loop inside a class
    FOR_ALL_CLASS(i, 0, nx_,
                  j, 0, ny_,
                  k, 0, nz_, {
                      
                  int idx = (i-1) * nx_ * ny_ + (j-1) * nz_ + (k-1);
                  arr3D_(i, j, k) = idx;
                  //printf("\nloop\n");
                  });
    Kokkos::fence();
    
    RUN_CLASS({
        printf("inside RUN_CLASS\n");
    });
    
}; // end member function



// functions called INSIDE a kokkos parallel loop
KOKKOS_INLINE_FUNCTION
void pass_by_ref(const FMatrixKokkos<int> &matrix);

KOKKOS_INLINE_FUNCTION
void pass_view_by_ref(const ViewFMatrixKokkos<int> &matrix);

KOKKOS_INLINE_FUNCTION
void pass_by_val(const FMatrixKokkos<int> matrix);



// functions NOT called in a kokkos parallel loop
void pass_by_ref_two(const FMatrixKokkos<int> &matrix);

FMatrixKokkos <int> return_by_val();




// function objects ModelA and ModelB with parallel loops inside
class ModelA{

private:

    // initial variables
    
public:

    // default constructor
    ModelA (){};
    
    // overload constructor to set initial values
    
    // overload()
    void operator()(FMatrixKokkos<int> &matrix){
        
	printf("inside ModelA \n");
	
        int loc_sum;
        int val = 0;
    
        // NOTE: if private vars are accessed, requires REDUCE_SUM_CLASS
        // do summation in parallel on GPU
        REDUCE_SUM_CLASS(k, 1, 6,
                         j, 1, 5,
                         i, 1, 4,
                         loc_sum, {
                         loc_sum += matrix(i,j,k,1);
                        }, val);
		  
	Kokkos::fence();

        printf(" val = %i \n", val);
	 
     } // end overload
    
}; // end function object ModelA



class ModelB{

private:

    // initial variables
    
public:

    // default constructor
    ModelB (){};
    
    // overload constructor to set initial values
    
    // overload()
    void operator()(FMatrixKokkos<int> &matrix){
        
	printf("inside ModelB \n");
	
        int loc_sum;
        int val = 0;
	
        // NOTE: if private vars are accessed, requires REDUCE_SUM_CLASS
        // do summation in parallel on GPU
        REDUCE_SUM_CLASS(k, 1, 6,
                         j, 1, 5,
                         i, 1, 4,
                         loc_sum, {
                         loc_sum += matrix(i,j,k,1);
                        }, val);
		  
	Kokkos::fence();

        printf(" val = %i \n", val);
	 
     } // end overload
    
}; // end function object ModelB
using models = std::variant<ModelA, ModelB>;



// function objects called inside a parallel loop

class MethodA{

private:

    // initial variables
    
public:

    // default constructor
    MethodA (){};
    
    // overload constructor to set initial values
    
    // overload()
    KOKKOS_INLINE_FUNCTION
    void operator()(const FMatrixKokkos<int> &matrix) const{
        
	printf("inside MethodA \n");
	
	int idx = matrix(1,1,1,1);
	
	matrix(1,1,1,1) = idx;  // do something pointless

    } // end overload
    
}; // end function object MethodA

class MethodB{

private:

    // initial variables
    
public:

    // default constructor
    MethodB (){};
    
    // overload constructor to set initial values
    
    // overload()
    KOKKOS_INLINE_FUNCTION
    void operator()(const FMatrixKokkos<int> &matrix) const{
        
	printf("inside MethodB \n");
	
	int idx = matrix(1,1,1,1);
	
	matrix(1,1,1,1) = idx;  // do something pointless

    } // end overload
    
}; // end function object MethodB


using methods = std::variant<MethodA, MethodB>;


/*
using my_variant = std::variant<std::monostate, A, B>;

void foo(my_variant &v) {
    switch (v.index()) {

    case 0: break; // do nothing because the type is std::monostate

    case 1: {
        doSomethingWith(std::get<A>(v));
        break;
    }

    case 2: {
        doSomethingElseWith(std::get<B>(v));
        break;
    }

    }
}
*/


template <typename F1, typename F2>
KOKKOS_INLINE_FUNCTION
void run_methods(const F1, const F2, methods);



// enum
namespace choices
{
    enum myChoice 
    { 
        METHOD_A = 1,
        METHOD_B = 2, 
        METHOD_C = 3  
    };
}



// function pointer
template <typename T>
struct method_ptrs{
   void (*fcn_ptr)(const T);
};

template <typename T>
KOKKOS_INLINE_FUNCTION
void sum(const T){
  printf("inside sum function\n");
};

template <typename T>
KOKKOS_INLINE_FUNCTION
void multiply(const T){
  printf("inside multiply function\n");
};



// struct that stores data inside
struct code_data_t{
   	double field_one[100];
 	int field_two[200];	
};


// a struct that stores MATAR dual arrays inside
struct cell_data_t{

    DCArrayKokkos <double> den;  
    DCArrayKokkos <double> pres;
    
    
    KOKKOS_INLINE_FUNCTION
    void initialize(const int i, const int j, const int k) const{
        den(i,j,k)  = 0.0;
	pres(i,j,k) = 0.0;
    };
    
};


// data in an exisiting framework that is managed
struct framework_data_t{
    
    // a 10 X 10 X 10 mesh
    int dim1 = 10;
    int dim2 = 10;
    int dim3 = 10;
    
    double data1[1000];  // notional data, could be dynammically allocated
    double data2[1000];  // ...
    
};

// view of data in an exisiting framework
struct framework_matar_t{

    DViewCArrayKokkos <double> data1;  // Views of the notional data on CPU and GPU
    DViewCArrayKokkos <double> data2;  // ...
    
};




//=============================================================
//
// Main function
//
int main(int argc, char *argv[]) {


    Kokkos::initialize(argc, argv);
    {   

        // -----------------------
        // parameters for examples
        // -----------------------
        u_int size_i, size_j, size_k, size_l;
        size_i = 3; size_j = 4; size_k = 5; size_l = 6;
        
        policy1D Arr_policy_1d = policy1D(0, size_i);
        policy2D Arr_policy_2d = policy2D({0, 0}, {size_i, size_j});
        policy3D Arr_policy_3d = policy3D({0, 0, 0}, {size_i, size_j, size_k});
        policy4D Arr_policy_4d = policy4D({0, 0, 0, 0}, {size_i, size_j, size_k, size_l});
        
        policy1D Mtx_policy_1d = policy1D(1, size_i+1);
        policy2D Mtx_policy_2d = policy2D({1, 1}, {size_i+1, size_j+1});
        policy3D Mtx_policy_3d = policy3D({1, 1, 1}, {size_i+1, size_j+1, size_k+1});
        policy4D Mtx_policy_4d = policy4D({1, 1, 1, 1}, {size_i+1, size_j+1, size_k+1, size_l+1});
        
        
        // -----------------------
        // CArray
        // -----------------------
        
        printf("\n1D CArray\n");
        auto cak1D = CArrayKokkos <int> (size_i);
        
        // a parallel 1D loop
        Kokkos::parallel_for("1DCArray", Arr_policy_1d, KOKKOS_LAMBDA(const int i) {
            cak1D(i) = i;
            //printf("%d) %d\n", i, cak1D(i));
        });
        Kokkos::fence();
        
        // the marco for a parallel 1D loop
        FOR_ALL(i, 0, size_i, {
            cak1D(i) = i;
        });
        
        Kokkos::fence();
        
        // -----------------------
        // FArray
        // -----------------------
        
        printf("\n2D FArray\n");
        auto fak2D = FArrayKokkos <int> (size_i, size_j);
        Kokkos::parallel_for("2DFArray", Arr_policy_2d, KOKKOS_LAMBDA(const int i, const int j) {
            int idx = j * size_i + i;
            fak2D(i, j) = idx;
            //printf("%d) %d\n", idx, fak2D(i, j));
        });
        Kokkos::fence();
        
        // the marco for a parallel 2D nested loop
        FOR_ALL(i, 0, size_i,
                j, 0, size_j,
                {
                int idx = j * size_i + i;
                fak2D(i, j) = idx;
                //printf("%d) %d\n", idx, fak2D(i, j));
                });
        Kokkos::fence();
        
	
        // -----------------------
        // CMatrix
        // -----------------------
        
        printf("\n3D CMatrix\n");
        auto cmk3D = CMatrixKokkos <int> (size_i, size_j, size_k);
	printf("made 3D CMatrix\n");
        printf("made CMATARkokkos\n");
        Kokkos::parallel_for("3DCMatrix", Mtx_policy_3d, KOKKOS_LAMBDA(const int i, const int j, const int k) {
            int idx = (i-1) * size_j * size_k + (j-1) * size_k + (k-1);
            cmk3D(i, j, k) = idx;
            printf("%d) %d\n", i, cmk3D(i, j, k));
        });
        Kokkos::fence();
        
        
        // the marco for a parallel 3D nested loop
        auto cmk3D_two = CMatrixKokkos <int> (size_i, size_j, size_k);
        FOR_ALL(i, 1, size_i+1,
                j, 1, size_j+1,
                k, 1, size_k+1,
                {
                int idx = (i-1) * size_j * size_k + (j-1) * size_k + (k-1);
                cmk3D_two(i, j, k) = idx;
                
                //printf("index %d) CMatrix = %d and %d\n", idx, cmk3D_two(i, j, k), cmk3D(i, j, k));
                });
        Kokkos::fence();
	
        
        // -----------------------
        // FMatrix
        // -----------------------
        
        printf("\n4D FMatrix\n");
        auto fmk4D = FMatrixKokkos <int> (size_i, size_j, size_k, size_l);
	
        Kokkos::parallel_for("4DFMatrix", Mtx_policy_4d, KOKKOS_LAMBDA(const int i, const int j, const int k, const int l) {
            int idx = (l-1) * size_i * size_j * size_k + (k-1) * size_i * size_j + (j-1) * size_i + (i-1);
            fmk4D(i, j, k, l) = idx;
        });
        Kokkos::fence();
	
	
	
	// -- functions exectuted on device inside a parallel for loop ---
	// A parallel loop
	FOR_ALL(i,0,1,{
	    pass_by_ref(fmk4D); 
	    pass_by_val(fmk4D);
	});
	Kokkos::fence();
	

	// --- call a function that has kokkos parallel loops inside it ---
	
	// get a FMatrix from a function
	fmk4D = return_by_val();
	
	// verify the values are correct
	FOR_ALL(i,0,1,{
	    pass_by_ref(fmk4D);
	});
	Kokkos::fence();    

	
	// call a function that has kokkos parallel loops inside it
	pass_by_ref_two(fmk4D);
	
	
	
	// -----------------------
        // ViewFMatrix
        // -----------------------
        
        printf("\n3D ViewFMatrix\n");
	
	// slice off the last dimension of FMatrix at L=1
	FOR_ALL(L,1,2,{
	    // MATAR views by default are on the device
            auto viewfmk3D = ViewFMatrixKokkos <int> (&fmk4D(1,1,1,L),size_i, size_j, size_k);
	    
	    // pass this view to a function
	    pass_view_by_ref(viewfmk3D);
	});
	Kokkos::fence();
	
	
	
	// -----------------------
        // functors
        // -----------------------
	printf("\nfunctors\n");
	ModelA model_a;
	model_a(fmk4D);
	
	// --
	MethodA method_a;
	FOR_ALL(i,1,2,{
	    method_a(fmk4D);
	});
	Kokkos::fence();  
	
	
	// -----------------------
        // std::variant access
        // -----------------------
	
	printf("\nstd::variant with functors\n");
		
        models my_model = ModelA(); // set model type
	
	size_t idx = my_model.index();
	printf("index of model in variant is = %lu \n",idx);

	// find and execute the model selected
	std::visit(overloaded {
            [&fmk4D](ModelA model) { 
	        printf("ModelA is being executed\n"); 
		
		model(fmk4D);

	    },
            [&fmk4D](ModelB model) { 
	        printf("ModelB is being executed\n"); 
		model(fmk4D);
	    }
        }, my_model);

	
	
	printf("\nCArray of std::variants with functors\n");
	// MATAR CArray of std::variants
	CArray <models> mat_models(3);
	mat_models(0) = ModelB(); // material 0 physics model
	mat_models(1) = ModelA(); // material 1 physics model
	mat_models(2) = ModelB(); // material 2 physics model
	
	idx = mat_models(0).index();
	printf("index of model in variant is = %lu \n",idx);
	
	for (int mat_id=0; mat_id<3; mat_id++){
	// find and execute the model selected
	std::visit(overloaded {
            [&fmk4D](ModelA model) { 
	        printf("ModelA is being executed\n"); 
		model(fmk4D);
	    },
            [&fmk4D](ModelB model) { 
	        printf("ModelB is being executed\n"); 
		model(fmk4D);
	    }
        }, mat_models(mat_id));
	} // end of loop over materials
	
	
	CArray <methods> mat_methods(3);
	mat_methods(0) = MethodB(); // material 0 numerical method
	mat_methods(1) = MethodA(); // material 1 numerical method
	mat_methods(2) = MethodB(); // material 2 numerical method

        // material centric approach
	for (int mat_id=0; mat_id<3; mat_id++){
	
	    // find and execute the model selected
	    std::visit(overloaded {
        	[&fmk4D](MethodA method) { 
	            printf("ModelA is being executed\n");

		    // e.g., loop over the cells in the mesh in parallel
		    FOR_ALL(i,1,2,{
	        	method(fmk4D);
	            });
		    Kokkos::fence();  

		},
        	[&fmk4D](MethodB method) { 
	            printf("ModelB is being executed\n"); 

		    // e.g., loop over the cells in the mesh in parallel
		    FOR_ALL(i,1,2,{
	        	method(fmk4D);
	            });
		    Kokkos::fence();  
		}
            }, mat_methods(mat_id));
	    
	} // end of loop over materials


        
	
	
	
    	// -----------------------
    	// DualView types
    	// -----------------------
	
	printf("\nDual views\n");
	
	code_data_t my_code_data;  // struct with arrays of data 

	
	// create a dual view of the data held inside my_code_data struct
	auto field_one = DViewCArrayKokkos <double> (&my_code_data.field_one[0], 100);
	auto field_two = DViewCArrayKokkos <int> (&my_code_data.field_two[0], 200);
	
	printf("modifying the dual view fields on the device\n");
	
	// modify the values in field one on the device 
    	FOR_ALL(i,0,100,{
	     field_one(i) = 12.345;
	});
	Kokkos::fence();
	field_one.update_host();  // copy data from devise to the host 
	
	printf("dual view of field_one = %f, struct field_one = %f \n", field_one.host(0), my_code_data.field_one[0]);
	
	// modify the values in field two on the device 
    	FOR_ALL(i,0,200,{
	     field_two(i) = 3;
	});
	Kokkos::fence();
	field_two.update_host();  // copy data from devise to the host 
	
	printf("dual view of field_two = %i, struct field_two = %i \n", field_two.host(0), my_code_data.field_two[0]);
	
	printf("modifying struct field_one = 314.5 \n");
	for (int i=0; i<100; i++){
	    my_code_data.field_one[i] = 314.15;
	} // end for loop
	printf("dual view of field_one = %f, struct field_one = %f \n", field_one.host(0), my_code_data.field_one[0]);
	


	// -----------------------
    // Dual Array types in an object
    // -----------------------
	
	printf("\nDual types inside struct\n");
	
	// struct with MATAR arrays of data 
	cell_data_t cell_data;  // allocate the data sizes: 10X10x10 mesh
        
	printf("allocate dual type sizes held in struct\n");
	cell_data.den  = DCArrayKokkos <double> (10,10,10);
	cell_data.pres = DCArrayKokkos <double> (10,10,10);

	
	// set the values inside the cell_data struct on the device
	
	printf("setting the dual type values and calling initialize functions \n");
        FOR_ALL(i, 0, 10,
                j, 0, 10,
                k, 0, 10,
                {
		   cell_data.initialize(i,j,k);
		   
                   cell_data.den(i,j,k) = 3.14159;
		   cell_data.pres(i,j,k) = 1.0;
                });
        Kokkos::fence();
	
	// update the host side
	cell_data.den.update_host();
	cell_data.pres.update_host();
	printf("The host values of the dual CArrays in the struct = %f and %f \n", cell_data.den.host(0,0,0), cell_data.pres.host(0,0,0));
	
	
	
	printf("\nDualView types inside struct\n");
	framework_data_t  framework_data;  // data is allocated by some framework across CPUs.
	
	// use MATAR to get the data onto the device e.g., GPU and make multiD views of the data
	framework_matar_t mtr_data;
	
	
	// get the mesh dims from the framework struct
	int mesh_dim1 = framework_data.dim1;
	int mesh_dim2 = framework_data.dim2;
	int mesh_dim3 = framework_data.dim3;
	
	printf("allocate data from the framework on the device\n");
	mtr_data.data1  = DViewCArrayKokkos <double> (&framework_data.data1[0], 
	    					      mesh_dim1, 
						      mesh_dim2,
						      mesh_dim3);
	mtr_data.data2  = DViewCArrayKokkos <double> (&framework_data.data2[0], 
	    					      mesh_dim1, 
						      mesh_dim2,
						      mesh_dim3);

		   
	printf("setting the dual type values\n");
	// set the framework values inside the struct on the device
        FOR_ALL(i, 0, mesh_dim1,
                j, 0, mesh_dim2,
                k, 0, mesh_dim3,
                {
                   mtr_data.data1(i,j,k) = 5.6;
		   mtr_data.data2(i,j,k) = 9.2;
                });
        Kokkos::fence();
	
	// update the host side
	mtr_data.data1.update_host();
	mtr_data.data2.update_host();
	printf("The 1st values of framework struct arrays = %f and %f \n", framework_data.data1[0], framework_data.data2[0]);
	// note how MATAR modified the data in the framework on the device
	
	// The dualView type also gives a view of the 1D framework data on the host side
	printf("The views of 1st host values of framework data  = %f and %f \n", mtr_data.data1.host(0,0,0), mtr_data.data2.host(0,0,0));
	
	framework_data.data1[0] = 77.77;
	framework_data.data1[1] = 88.88;
	mtr_data.data1.update_device();
	RUN({
	    printf("value on device after update = %f, %f", mtr_data.data1(0,0,0), mtr_data.data1(0,0,1));
	});
	Kokkos::fence();
	
	printf("\n");	
	
	
	
        
        // -----------------------
        // DynamicRaggedRightArray
        // -----------------------
        
        printf("\nDynamic Ragged Right Array\n");
        DynamicRaggedRightArrayKokkos <int> drrak;
        drrak = DynamicRaggedRightArrayKokkos <int> (size_i, size_j);
        
        Kokkos::parallel_for("DRRAKTest", size_i, KOKKOS_LAMBDA(const int i) {
            for (int j = 0; j < (i % size_j) + 1; j++) {
                drrak.stride(i)++;
                drrak(i,j) = j;
                //printf("(%i) stride is %d\n", i, j);
            }
        });
        Kokkos::fence();
	
        printf("\ntesting macro FOR_ALL\n");
	
        // testing MATAR FOR_ALL loop
        DynamicRaggedRightArrayKokkos <int> my_dyn_ragged(size_i, size_j);
        FOR_ALL(i, 0, size_i, {
            for (int j = 0; j <= (i % size_j); j++) {
                my_dyn_ragged.stride(i)++;
                my_dyn_ragged(i,j) = j;
		printf(" dyn_ragged_right error = %i \n", my_dyn_ragged(i,j)-drrak(i,j));
            }// end for
        });// end parallel for
        Kokkos::fence();
        
        // -----------------------
        // RaggedRightArray
        // -----------------------
        printf("\nRagged Right Array\n");
        // testing ragged initialized with CArrayKokkos for strides
        CArrayKokkos <size_t> some_strides(4);
        
        // create a lower-triangular array
        RUN({
            some_strides(0) = 1;
            some_strides(1) = 2;
            some_strides(2) = 3;
            some_strides(3) = 4;
        });
        
        RaggedRightArrayKokkos <int> lower_tri(some_strides);
        
        
        
        
        // -----------------------
        // CArray view
        // -----------------------
        
        printf("\nView CArray\n");
        std::array<int, 9> A1d;
        for (int init = 0; init < 9; init++) {
            A1d[init] = init+1;
        }
        policy2D CAKPpol = policy2D({0,0}, {3, 3});
        DViewCArrayKokkos <int> cakp;
        cakp = DViewCArrayKokkos <int> (&A1d[0], 3, 3);
        Kokkos::parallel_for("CAKPTest", CAKPpol, KOKKOS_LAMBDA(const int i, const int j) {
            //printf("%d) %d\n", i * 3 + j, cakp(i, j));
        });
        Kokkos::fence();
        
        // -----------------------
        // CArray inside a class
        // -----------------------
        
        printf("\nCArray in a class\n");
        Data my_data(size_i, size_j, size_k);
        my_data.some_fcn();
	
	
	
	
	printf("\nENUM\n");
	
	// simple enum example:
	//    choices::myChoice enumVar;
	//    enumVar = choices::METHOD_A; // setting the method
	
	// declare methods
	MethodA my_method_a;
	MethodB my_method_b;

		
	printf("CArrayKokkos of enums\n");
	auto time_1 = std::chrono::high_resolution_clock::now();
	CArrayKokkos <choices::myChoice> my_choices(2);
	
	// set the method on the GPU
	RUN({
	    my_choices(0) = choices::METHOD_A;
	    my_choices(1) = choices::METHOD_B;
	});
	Kokkos::fence();  
	
	
	
 	// e.g., loop over in parallel
	FOR_ALL(i,1,2,{
	    printf("selecting method\n");
	    
	    switch (my_choices(i))
            {
        	case choices::METHOD_A:
        	{
        	    // do stuff
		    printf("using method_A\n");
		    my_method_a(fmk4D);
        	    break;
        	}
		  
        	case choices::METHOD_B:
        	{
        	    // do stuff
		    printf("using method_B\n");
		    my_method_b(fmk4D);
        	    break;
        	}
		
        	default:
        	{
        	  // do nothing
        	}
            };  // end switch
	    
	    
	});
	Kokkos::fence();  
	
	auto time_2 = std::chrono::high_resolution_clock::now();
	
	
	std::cout << "Elapsed time in seconds: "
                  << std::chrono::duration_cast<std::chrono::microseconds>(time_2 - time_1).count()
                  << " microsec" << std::endl;
	
        
	
	printf("\nCArray of function pointers\n");
	
	//method_ptrs;
	CArrayKokkos < method_ptrs<FMatrixKokkos<int>> > Array_ptrs(2);
	
	
	// set the pointer on the device e.g., GPU
	RUN({
	    Array_ptrs(0).fcn_ptr = sum;
	    Array_ptrs(1).fcn_ptr = multiply;
	});
	Kokkos::fence();  
	
	// use the function
	RUN({
	    Array_ptrs(0).fcn_ptr(fmk4D);
	    Array_ptrs(1).fcn_ptr(fmk4D);
	});
	Kokkos::fence();
	
	
    } // end of kokkos scope
    
    Kokkos::finalize();

    printf("\nfinished\n\n");

    return 0;
}



// -----  Functions called INSIDE a kokkos parallel loop -----

KOKKOS_INLINE_FUNCTION
void pass_by_ref(const FMatrixKokkos<int> &matrix){
  printf("inside pass_by_ref function,");


  int val = 0;
  for (int k=1; k<=5; k++){
    for (int j=1; j<=4; j++){
      for (int i=1; i<=3; i++){
		
        val += matrix(i,j,k,1);
		
      } // end for i
    } // end for j
  } // end for k

  printf(" val = %i \n", val);

}


KOKKOS_INLINE_FUNCTION
void pass_by_val(const FMatrixKokkos<int> matrix){
  printf("inside pass_by_val function,");

  int val = 0;
  
  // do summation in serial
  for (int k=1; k<=5; k++){
    for (int j=1; j<=4; j++){
      for (int i=1; i<=3; i++){
		
        val += matrix(i,j,k,1);
		
      } // end for i
    } // end for j
  } // end for k

  printf(" val = %i, \n", val);

}


KOKKOS_INLINE_FUNCTION
void pass_view_by_ref(const ViewFMatrixKokkos<int> &matrix){

  // remember that MATAR views are always on the device
  
  printf("inside pass_view_by_ref function,");

  int val = 0;
  
  // do summation in serial
  for (int k=1; k<=5; k++){
    for (int j=1; j<=4; j++){
      for (int i=1; i<=3; i++){
		
        val += matrix(i,j,k);
		
      } // end for i
    } // end for j
  } // end for k

  printf(" val = %i, \n", val);

} // end function


template <typename F1, typename F2>
KOKKOS_INLINE_FUNCTION
void run_methods(const F1 &lambda_fcn1, const F2 &lambda_fcn2, methods &v) {
    switch (v.index()) {

        case 0: {
            lambda_fcn1(std::get<MethodA>(v));
            break;
        } // end case 0

        case 1: {
            lambda_fcn2(std::get<MethodB>(v));
            break;
        } // end case 1

    } // end case 
}; // end of function






// -----  Functions NOT called in a kokkos parallel loop -----

void pass_by_ref_two(const FMatrixKokkos<int> &matrix){
  printf("inside pass_by_ref_two function (parallel loops),");

    int loc_sum;
    int val = 0;
    
    // do summation in parallel on GPU
    REDUCE_SUM(k, 1, 6,
               j, 1, 5,
               i, 1, 4,
               loc_sum, {
                   loc_sum += matrix(i,j,k,1);
               }, val);

  printf(" val = %i \n", val);

}


FMatrixKokkos <int> return_by_val(){

        printf("inside return_by_val \n");

        // -----------------------
        // parameters for examples
        // -----------------------
        u_int size_i, size_j, size_k, size_l;
        size_i = 3; size_j = 4; size_k = 5; size_l = 6;
        
        policy4D Mtx_policy_4d = policy4D({1, 1, 1, 1}, {size_i+1, size_j+1, size_k+1, size_l+1});

        FMatrixKokkos <int> fmk4D_local(size_i, size_j, size_k, size_l);
	
        Kokkos::parallel_for("4DFMatrix", Mtx_policy_4d, KOKKOS_LAMBDA(const int i, const int j, const int k, const int l) {
            int idx = (l-1) * size_i * size_j * size_k + (k-1) * size_i * size_j + (j-1) * size_i + (i-1);
            fmk4D_local(i, j, k, l) = idx;
        });
        Kokkos::fence();
	
	return fmk4D_local;
}





