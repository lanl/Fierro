# SGH-RZ Solver

The SGH-RZ solver is an axisymmetric Lagrangian hydrodynamic method designed for 2D axisymmetric coordinates $(z, r, \varphi)$ (Figure~\ref{fig:rz-coord-system}). It conserves mass and total energy and preserves symmetry on 1D radial flows using equal angle polar meshes.

## Governing Equations

The solver addresses the conservation of mass, momentum, and specific internal energy in cylindrical coordinates:

1.  **Mass Conservation:**
    The density for each material is calculated using strong mass conservation.
    $$
    \rho^m_h = \frac{m^m_h}{V^m_h}
    $$
    where $V^m_h = \theta^m_h \beta_h^m V^{3D}_h$.

2.  **Momentum Evolution:**
    The Petrov-Galerkin approach is used with a specific test function $\eta_q = \phi_q \frac{r_q}{r}$ to preserve symmetry:
    $$
    \sum \limits_{h \ni q} \sum \limits_{m \in h} \int \limits_{V_h} \eta_q \rho^m_h \frac{d {\bf v}_h}{dt} \beta_h^m \theta_h^m dV = \sum \limits_{h \ni q} \sum \limits_{m \in h} \int \limits_{V_h} \eta_q \nabla \cdot (\boldsymbol{\sigma}_h^m + {\bf Q}_h^m) \beta_h^m \theta_h^m dV
    $$

3.  **Specific Internal Energy Evolution:**
    The specific internal energy evolution equation guarantees total energy conservation (compatible discretization). The change in specific internal energy for an element is given by:
    $$
    {e}_h^{m,\, n+1} = {e}_h^{m,\, n} - \frac{\Delta t}{m^m_h} \sum \limits_{p \in h} \left( { {\bf F}^{m, \,n+1/2}_{hp}}\cdot \frac{1}{2}\left({\bf v}_p^{n+1} + {\bf v}_p^{n} \right) \right)
    $$

## Numerical Methods

### Geometry
Position and velocity fields are defined in terms of Lagrangian basis functions in 2D:
$$
{\bf x}_h({\boldsymbol  \xi},t) = \sum \limits_{p \in h} {\phi}_p \left( {\boldsymbol  \xi} \right) \cdot {\bf x}_p \left( t \right)
$$
$$
{\bf v}_h({\boldsymbol  \xi},t) = \sum \limits_{p \in h} {\phi}_p \left( {\boldsymbol  \xi} \right) \cdot {\bf v}_p \left( t \right)
$$

![RZ Coordinate System](../../../../../docs/figs/meshRZcoordinates.png)
![RZ Element Mapping](../../../../../docs/figs/LinearOrderRZMappings.png)

### Time Integration
The discrete change in velocity is calculated using a two-step Runge-Kutta time integration method:
1.  Calculate intermediate velocity ${\bf v}_p^{n+1/2}$.
2.  Calculate final velocity ${\bf v}_p^{n+1}$.

### Artificial Viscosity
A tensoral dissipation term ${\bf Q}^m_h$ is included for stability on shock problems, calculated using a multi-directional approximate Riemann solver (MARS).

### Stress Components
![Stress Components](../../../../../docs/figs/stressRZelement.png)
The stress tensor includes radial ($r$), axial ($z$), and azimuthal ($\varphi$) components.

## Implementation Details
-   **Element Type:** Single quadrature point quadrilateral element in 2D representing a revolved volume.
-   **Mass Lumping:** Conservative and consistent partitioning of element area to corners.
-   **Material Handling:** Supports arbitrary number of materials per element with volume fractions $\theta^m_h$ and material fractions $\beta_h^m$.
