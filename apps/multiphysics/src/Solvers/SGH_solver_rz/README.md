# SGH-RZ Solver

The SGH-RZ solver is an axisymmetric Lagrangian hydrodynamic method designed for 2D axisymmetric coordinates $(z, r, \varphi)$. It conserves mass and total energy and preserves symmetry on 1D radial flows using equal angle polar meshes.

![RZ Coordinate System](../../../../../docs/figs/meshRZcoordinates.png)

## Governing Equations

Even though the simulation is performed using a 2D mesh, the simulation corresponds to a 3D problem that is axisymmetric. The axisymmetric approximation applies to cases with:

$$
\frac{\partial \rho^{3D}} {\partial \varphi} = 0, \quad \frac{\partial {\mathbf v}^{3D}} {\partial \varphi} = 0, \quad \text{and} \quad \frac{\partial e^{3D}} {\partial \varphi} = 0
$$

The superscript 3D will be used to differentiate a variable from one that is purely planar. The axisymmetric approximation reduces a 3D system $(z,r,\varphi)$ to a 2D solve $(z,r)$.


## Numerical Methods

The discussion that follows here will focus on a novel numerical formulation to solve a subset of the governing equations presented in this section for an arbitrary number of materials in an element.

### Mass Conservation
The density for each material is calculated using strong mass conservation.

$$
\rho^m_h = \frac{m^m_h}{V^m_h}
$$

where:
- $\rho^m_h$ is the density of material $m$ in element $h$.
- $m^m_h$ is the mass of material $m$ in element $h$.
- $V^m_h$ is the volume of material $m$ in element $h$.

The material volume is defined as:

$$
V^m_h = \theta^m_h \beta_h^m V^{3D}_h
$$

where:
- $\theta^m_h$ is the volume fraction of material $m$ in element $h$.
- $\beta_h^m$ is the material fraction of material $m$ in element $h$.
- $V^{3D}_h$ is the total volume of element $h$.

### Momentum Evolution
The Petrov-Galerkin approach is used with a specific test function $\eta_q = \phi_q \frac{r_q}{r}$ to preserve symmetry:

$$
\sum \limits_{h \ni q} \sum \limits_{m \in h} \int \limits_{V_h} \eta_q \rho^m_h \frac{d \mathbf{v}_h}{dt} \beta_h^m \theta_h^m dV = \sum \limits_{h \ni q} \sum \limits_{m \in h} \int \limits_{V_h} \eta_q \nabla \cdot (\boldsymbol{\sigma}_h^m + \mathbf{Q}_h^m) \beta_h^m \theta_h^m dV
$$

where:
- $q$ represents the node index.
- $h$ represents the element index.
- $\eta_q$ is the test function associated with node $q$.
- $\phi_q$ is the basis function associated with node $q$.
- $r_q$ is the radial position of node $q$.
- $r$ is the radial position within the element.
- $\mathbf{v}_h$ is the velocity field within element $h$.
- $\boldsymbol{\sigma}_h^m$ is the Cauchy stress tensor for material $m$ in element $h$.
- $\mathbf{Q}_h^m$ is the artificial viscosity tensor for material $m$ in element $h$.

### Specific Internal Energy Evolution
The specific internal energy evolution equation guarantees total energy conservation (compatible discretization). The change in specific internal energy for an element is given by:

$$
e^{m, n+1/2}_h =e^{m, \, n}_h - \frac{1}{2}\frac{\Delta t}{m^{m, \, 3D}_h}\sum \limits_{hp \in h}  {\bf F}^m_{hp}r_p \cdot \frac{1}{2}\left( {\bf v}_p^{n+1/2} + {\bf v}_p^{n} \right).
$$

$$
e^{m,  n+1}_h =e^{m, n+1/2}_h - \frac{\Delta t}{m^{m, \, 3D}_h}\sum \limits_{hp \in h}  {\bf F}^m_{hp}r_p \cdot 
    \frac{1}{2}\left( {\bf v}_p^{n+1} + {\bf v}_p^{n} \right).
$$

where:
- ${e}_h^{m, n+1}$ is the specific internal energy of material $m$ in element $h$ at time step $n+1$.
- ${e}_h^{m, n}$ is the specific internal energy at time step $n$.
- $\Delta t$ is the time step size.
- $p$ represents the node index of element $h$.
- $\mathbf{F}^{m, \,n+1/2}_{hp}$ is the corner force exerted by material $m$ in element $h$ on node $p$ at the half time step.
- $\mathbf{v}_p^{n+1}$ and $\mathbf{v}_p^{n}$ are the velocities of node $p$ at time steps $n+1$ and $n$, respectively.

### Geometry
Position and velocity fields are defined in terms of Lagrangian basis functions in 2D:

$$
\mathbf{x}_h({\boldsymbol \xi},t) = \sum \limits_{p \in h} {\phi}_p \left( {\boldsymbol \xi} \right) \cdot \mathbf{x}_p \left( t \right)
$$

$$
\mathbf{v}_h({\boldsymbol \xi},t) = \sum \limits_{p \in h} {\phi}_p \left( {\boldsymbol \xi} \right) \cdot \mathbf{v}_p \left( t \right)
$$

where:
- $\mathbf{x}_h({\boldsymbol \xi},t)$ is the position field within element $h$.
- $\mathbf{v}_h({\boldsymbol \xi},t)$ is the velocity field within element $h$.
- ${\boldsymbol \xi}$ represents the reference coordinates.
- ${\phi}_p$ is the basis function associated with node $p$.
- $\mathbf{x}_p(t)$ is the position of node $p$ at time $t$.
- $\mathbf{v}_p(t)$ is the velocity of node $p$ at time $t$.

![RZ Element Mapping](../../../../../docs/figs/LinearOrderRZMappings.png)

### Time Integration
The discrete change in velocity is calculated using a two-step Runge-Kutta time integration method:
1.  Calculate intermediate velocity $\mathbf{v}_p^{n+1/2}$.
2.  Calculate final velocity $\mathbf{v}_p^{n+1}$.

### Artificial Viscosity
A tensoral dissipation term $\mathbf{Q}^m_h$ is included for stability on shock problems, calculated using a multi-directional approximate Riemann solver (MARS).

## Implementation Details
-   **Element Type:** Single quadrature point quadrilateral element in 2D representing a revolved volume.
-   **Mass Lumping:** Conservative and consistent partitioning of element area to corners.
-   **Material Handling:** Supports arbitrary number of materials per element with volume fractions $\theta^m_h$ and material fractions $\beta_h^m$.
