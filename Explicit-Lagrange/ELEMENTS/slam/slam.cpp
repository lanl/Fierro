#include <iostream>  // std::cout etc.
#include <math.h> 


#include "utilities.h"
#include "matar.h"
#include "slam.h"


using namespace utils;



////////////////////////////
// Function Definitions   //
////////////////////////////


/* ------------------------- */
/* LU decomposition function */
/* ------------------------- */
int LU_decompos(
    ViewCArray <real_t> &source_mat,
    ViewCArray <int> &indx,           // permutations
    int &parity,                        // parity (+1 0r -1)
    const int n) {                      // size of matrix


    int i, imax, j, k; // Indexing 

    real_t big, sum, temp;// useful 

    real_t vv[n];       // temp arrary for solver

    parity = 1;
    
    /* search for the largest element in each row; save the scaling in the 
    temporary array vv and return zero if the matrix is singular */
    for(i = 0; i < n; i++) {
        
        big = 0.;
        for(j = 0; j < n; j++) if((temp=fabs(source_mat(i,j))) > big) big=temp;
        
        if(big == 0.) return(0);
        
        vv[i] = big;
    }

    /* the main loop for the Crout's algorithm */
    for(j = 0; j < n; j++) {
        
        /* this is the part a) of the algorithm except for i==j */
        for(i=0;i<j;i++) {
            
            sum=source_mat(i,j);
            
            for(k=0;k<i;k++) sum -= source_mat(i,k)*source_mat(k,j);

            source_mat(i,j) = sum;
        }
    
        /* initialize for the search for the largest pivot element */
        big = 0.;
        imax = j;
        
        /* this is the part a) for i==j and part b) for i>j + pivot search */
        for(i = j; i < n; i++) {
            
            sum = source_mat(i,j);
            
            for(k=0; k<j; k++) sum -= source_mat(i,k)*source_mat(k,j);
            
            source_mat(i,j) = sum;
            
            /* is the figure of merit for the pivot better than the best so far? */
            if((temp = vv[i]*fabs(sum)) >= big) {big = temp; imax = i;}
        }

        /* interchange rows, if needed, change parity and the scale factor */
        if(imax != j) {
            
            for(k = 0; k < n; k++){
                temp = source_mat(imax,k);
                source_mat(imax,k) = source_mat(j,k);
                source_mat(j,k) = temp;
            }
            
            parity = -(parity);
            vv[imax] = vv[j];
        }
        
        /* store the index */
        indx(j) = imax;
        /* if the pivot element is zero, the matrix is singular but for some 
        applications a tiny number is desirable instead */
        
        if(source_mat(j,j) == 0.) source_mat(j,j) = TINY;
        /* finally, divide by the pivot element */
        
        if(j<n-1) {
            
            temp=1./source_mat(j,j);
            for(i = j+1; i < n; i++) source_mat(i,j)*=temp;
        }
    }
    
    return(1);
}


/* ----------------------------- */
/* LU back substitution function */
/* ----------------------------- */

void LU_backsub(
    ViewCArray <real_t> &source_mat,  // input matrix
    const ViewCArray <int> &indx,        // permutations
    ViewCArray <real_t> &b,      // Least squares coefficents
    const int n) {
    



        int i, j, ip, ii = -1;
        double sum;
        
        /* First step of backsubstitution; the only wrinkle is to unscramble 
        the permutation order. Note: the algorithm is optimized for a 
        possibility of large amount of zeroes in b */
        
        for(i = 0; i < n; i++) {
            
            ip = indx(i);
            // std::cout<<" ip = "<< ip <<std::endl;
            sum = b(ip);
            b(ip) = b(i);
            
            // source_mat(i,j)
            
            if(ii >= 0) for(j = ii; j<i; j++) sum -= source_mat(i,j)*b(j);
            else if(sum) ii=i;  /* a nonzero element encounted */
            
            b(i) = sum;
        }
        
        /* the second step */
        for(i=n-1; i>=0; i--) {
            
            sum = b(i);
            for(j=i+1; j<n; j++) sum-=source_mat(i,j)*b(j);
            
            b(i)=sum/source_mat(i,i);
        }

        
}


/* ------------------ */
/* LU invert function */
/* ------------------ */
void LU_invert(
    ViewCArray <real_t> &source_mat,  // input matrix
    ViewCArray <int> &indx,        // permutations
    ViewCArray <real_t> &inv_mat,     // inverse matrix
    ViewCArray <real_t> &col,         // tmp array
    int n) {


    for(int j = 0; j < n; j++){
        for(int i = 0; i < n; i++) col(i) = 0.0;
        
        col(j) = 1.0;
        LU_backsub(source_mat, indx, col, n);
        
        for(int i = 0; i < n; i++) inv_mat(i,j) = col(i);
    }
}


/* ----------------------- */
/* LU determinant function */
/* ----------------------- */
double LU_determ(
    ViewCArray <real_t> &source_mat,  // input matrix
    ViewCArray <int> &indx,        // permutations
    const int parity,                   // parity (+1 0r -1)
    const int n) {

    int j;
    real_t res = (real_t)(parity);
    
    for(j=0; j<n; j++) res *= source_mat(j,j);

    return(res);
}


