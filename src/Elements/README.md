# ELEMENTS

[![Linux](https://github.com/lanl/ELEMENTS/actions/workflows/Linux.yaml/badge.svg)](https://github.com/lanl/ELEMENTS/actions/workflows/Linux.yaml)
[![MacOS](https://github.com/lanl/ELEMENTS/actions/workflows/MacOS.yaml/badge.svg)](https://github.com/lanl/ELEMENTS/actions/workflows/MacOS.yaml)

## What is ELEMENTS?

The C++ **ELEMENTS** library is a collection of sub-libraries to support implementing a diverse range of numerical methods on low and high-order meshes.  The **ELEMENTS** library can be used for research and development of both continuous and discontinuous finite element methods, as well as, finite volume methods to solve a diverse range of partial differential equations. The **ELEMENTS** library includes the following sub-libraries:  **MATAR** contains the routines to support dense and sparse **mat**rices and **ar**rays, **SLAM** contains the interfaces to **s**olvers, **l**inear **a**lgebra, and **m**athematical routines or external packages (e.g., Trilinos),  **elements** contains the mathematical functions to support a large range of elements types including serendipity elements, **SWAGE** contains the routines and data-structures to support unstructured arbitrary-order 3D meshes that move or remain stationary, and **geometry** combines together **SWAGE** and **elements**.  The **ELEMENTS** libary is designed to support Lagrangian (mesh moves) solid dynamics and mechanics codes, Eulerian (mesh is stationary) fluid dynamics codes, and many other code applications.  

<p align="center"><img src="https://github.com/lanl/ELEMENTS/blob/develop-msu/docs/images/codeStructureELEMENTS.png" width="400">
<p align="center">Fig. Code structure layout
  
<p align="center"><img src="https://github.com/lanl/ELEMENTS/blob/develop-msu/docs/images/TaylorGreenVortex-t0.png" width="400"><img src="https://github.com/lanl/ELEMENTS/blob/develop-msu/docs/images/TaylorGreenVortex-tEnd.png" width="400">
<p align="center">Fig. A high-order 3D mesh deforming in the Taylor-Green vortex

## Getting started

To learn more about ELEMENTS and how to get started using it, please see the [ELEMENTS documentation](https://lanl.github.io/ELEMENTS/).

## How to cite

If you use the ELEMENTS library in your work, please cite the following in any pursuant research papers.

```
@article{MOORE2019100257,
  title = "{ELEMENTS: A high-order finite element library in C++}",
  journal = {SoftwareX},
  volume = {10},
  pages = {100257},
  year = {2019},
  issn = {2352-7110},
  doi = {https://doi.org/10.1016/j.softx.2019.100257},
  url = {https://www.sciencedirect.com/science/article/pii/S235271101930113X},
  author = {Jacob L. Moore and Nathaniel R. Morgan and Mark F. Horstemeyer},
  keywords = {Element Library, C++, High-order elements, Spectral elements, Serendipity elements}
}
```
