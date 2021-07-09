# Fierro

**FIERRO** (LANL code number C21030) is a modern C++ code intended to simulate quasi-static solid mechanics problems and transient, compressible material dynamic problems with Lagrangian methods, which have meshes with constant mass elements that move with the material.  **FIERRO** is designed to aid a) material model research that has historically been done using commercial implicit and explicit finite element codes, b) numerical methods research, and c) computer science research.  **FIERRO** supports user developed material models that adhere to several industry standard formats by using a C++ to Fortran interface to couple the model to the numerical solvers.  **FIERRO** is built on the **ELEMENTS** library that supports a diverse suite of element types, including high-order elements, and quadrature rules. The mesh class within the **ELEMENTS** library is designed for efficient calculations on unstructured meshes and to minimize memory usage.  **FIERRO** is designed to readily accommodate a range of numerical methods including continuous finite element, finite volume, and discontinuous Galerkin methods.  **FIERRO** is designed to support explicit and implicit time integration methods.  


## Computer implementation
**FIERRO** is implemented in C++ following modern programming practices.  **FIERRO** leverages the unique features of the **ELEMENTS** library, as such, this code serves as an example for solving a system of partial differential equations using the mesh class and geometric functions within the **ELEMENTS** library.  **FIERRO** registers state at material points within the element, registers polynomial fields in the element, and registers kinematic variables at element vertices.  The routines for the state are designed for highly efficient computations and to minimize memory usage.  The loops are written to aid fine-grained parallelization and to allow vectorization. **FIERRO** is a light-weight software application, and cleanly written following modern programming practices, so it useful for researching computer science based technologies for software performances portability.  

## Spatial discretization methods 
**FIERRO** has an established conservative Lagrangian finite element method for solving the governing equations (mass, momentum, and specific internal energy evolution equations) for compressible material dynamics using unstructured linear hexahedral elements combined with a multidirectional approximate Riemann solver for improved accuracy on smooth flows and stable solutions near velocity discontinuities and large gradients in a flow. **FIERRO** is designed for Lagrangian methods research and development so other types of numerical methods can be readily added to the code including Lagrangian finite volume and discontinuous Galerkin methods.  Likewise, high-order Lagrangian methods can be studied within the code because it is built upon the **ELEMENTS** library that supports high-order elements and high-order quadrature rules.  Numerical methods are being added to **FIERRO** to simulate quasi-static solid mechanics problems.

## Temporal discretization methods 
**FIERRO** supports a range of published multi-step time integration methods.  The code has an explicit multi-step Runge Kutta time integration method.  Implicit time integration methods can be implemented in **FIERRO**.

## Material models  
The classical ideal gas model is the only material model implemented in the code, and it is useful for verification tests of the software and simple gas dynamic simulations.  The code has C++ to Fortran interfaces that allow users to write their own material models and test them on quasi static problems or material dynamic applications.  The interfaces follow an industry standard format so that **FIERRO** can be used for model research and development that has historically been done with commercial implicit or explicit finite element codes. 














