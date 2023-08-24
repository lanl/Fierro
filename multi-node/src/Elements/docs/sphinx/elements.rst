.. _elements:

elements
========
A mesh is composed of non-overlapping elements, where the **elements** sub-library contains the mathematical functions to support a very broad range of element types including: 

* linear, quadratic, and cubic serendipity elements in 2D and 3D; 
* arbitrary-order tensor product elements in 2D and 3D;
* arbitrary-order spectral elements; and 
* a linear 4D element. 

The **elements** sub-library has functions to calculate quantities that are commonly used in finite element methods (both continuous and discontinuous) such as a basis function, gradient of a basis function, the Jacobian matrix, the inverse Jacobian matrix, the determinant of the Jacobian matrix, and a physical position inside the element, to name a few examples. 
The **elements** sub-library also supports both Gauss-Legendre and Gauss-Lobatto quadrature rules up to 8 quadrature points in each coordinate direction. 

elements
--------
A mesh is a collection of non-overlapping elements, and the mesh can be unstructured.  
Each reference element is defined in a reference coordinate system with a spatial mapping (e.g., an interpolation polynomial) to the physical coordinates.  
An element in the physical coordinates can have edges that are defined by straight lines (i.e., linear) or have edges defined by a high-order polynomial. 
The **elements** sub-library contains the mathematical functions for a suite of element types.    
The **geometry** sub-library combines the mesh data structures and mesh index spaces in **SWAGE** with the reference element defined in **elements**. 

geometry definition
-------------------
The position inside an element is given by a polynomial.  
**elements** contains a suite of spatial polynomials and also a space-time polynomial.  
The Jacobi matrix is equal to the gradient of this spatial polynomial.  
The volume of an element can be calculated using the determinate of the Jacobi matrix.  
The inverse of the Jacobi matrix is used to transform a Nabla operator from the physical coordinates to the reference coordinates.

index conventions
-----------------
The reference element contains nodes that are collocated with the Lobatto quadrature points.  
The local indexing of the quadrature points and nodes within an element follows an (i,j,k) index convention. 
An index map is supplied to convert the Serendipity local index convention to the local (i,j,k) index convention.

.. doxygennamespace:: elements
   :members:
