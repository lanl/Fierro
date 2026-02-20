# EVPFFT Micromechanical Solvers

Fierro provides a suite of Elasto-Viscoplastic Fast Fourier Transform (EVPFFT) micromechanical solvers designed for computationally efficient modeling of polycrystalline materials. These solvers leverage the spectral method to solve the governing equations of equilibrium and compatibility on a regular voxel grid, enabling the simulation of complex microstructures with high fidelity.

For detailed information, please refer to the paper: **Fierro, Enabling Multi-Material Multi-Physics**.

## Supported Solver Types

The following EVPFFT solver variants are available in Fierro, each tailored to specific modeling needs:

### 1. EVPFFT (Small Strain)
The standard EVPFFT implementation for small strain approximations.
*   **Description:** Solves the governing equations under the assumption of infinitesimal strains. It uses a Green's function method combined with an Augmented Lagrangian scheme to enforce compatibility and equilibrium.
*   **Use Cases:** Ideal for problems where deformations are small, such as elastic loading, early stages of plasticity, or fatigue analysis where geometric non-linearities are negligible. Efficient for crystal plasticity modeling of polycrystalline aggregates.

### 2. LS-EVPFFT (Large Strain)
An extension of the EVPFFT framework to finite strains (large deformations).
*   **Description:** Formulated in a Lagrangian frame (or updated Lagrangian), this solver accounts for geometric non-linearities, including grain rotation and texture evolution during large plastic deformations. It solves for the velocity gradient field rather than the displacement gradient.
*   **Use Cases:** Necessary for simulating metal forming processes (rolling, extrusion), high-velocity impacts, and any application involving significant plastic deformation where texture evolution and geometric changes are critical.

### 3. LSNP-EVPFFT (Large Strain Non-Periodic)
A variation of the Large Strain EVPFFT solver that relaxes the strict periodic boundary conditions.
*   **Description:** While standard FFT-based methods inherently assume periodic boundary conditions, LSNP-EVPFFT incorporates modifications to handle non-periodic conditions, often through techniques like buffer zones or modified Green's functions.
*   **Use Cases:** Suitable for simulations where periodic boundaries are not physical, such as modeling free surfaces, specific experimental setups, or problems with non-periodic loading conditions.

### 4. LS-EVPFFT-J2 (Large Strain J2 Plasticity)
A specialized version of the Large Strain solver using J2 (von Mises) plasticity.
*   **Description:** Instead of crystal plasticity models that track slip systems in individual grains, this solver uses a phenomenological J2 plasticity model (isotropic hardening).
*   **Use Cases:** Useful for comparing crystal plasticity results with continuum plasticity models, or for simulating materials where a J2 approximation is sufficient (e.g., isotropic metals, certain polymers) within a spectral framework.

## Green's Function Method

The core of the EVPFFT approach relies on transforming the differential equations of equilibrium into integral equations using Green's functions.

### Small Strain Formulation
The Cauchy stress $\boldsymbol{\sigma}(\mathbf{x})$ is decomposed into a polarization field $\boldsymbol{\phi}(\mathbf{x})$ relative to a homogeneous reference medium with stiffness $\mathbf{C}^0$:
$$
\boldsymbol{\sigma}(\mathbf{x}) = \mathbf{C}^{0}:\nabla\mathbf{u}(\mathbf{x}) + \boldsymbol{\phi}(\mathbf{x})
$$
where $\boldsymbol{\phi}(\mathbf{x}) = \boldsymbol{\sigma}(\mathbf{x}) - \mathbf{C}^{0}:\nabla\mathbf{u}(\mathbf{x})$.

Substituting this into the equilibrium equation $\nabla \cdot \boldsymbol{\sigma} = \mathbf{0}$ yields:
$$
\mathbf{C}^{0}:\nabla \nabla \mathbf{u}(\mathbf{x}) + \nabla \cdot \boldsymbol{\phi}(\mathbf{x}) = \mathbf{0}
$$

The solution for the displacement gradient $\nabla\mathbf{u}(\mathbf{x})$ is given by the convolution of the periodic Green's operator $\boldsymbol{\Gamma}$ and the polarization field:
$$
\nabla\mathbf{u}(\mathbf{x}) = \mathbf{E} - \boldsymbol{\Gamma}(\mathbf{x}) * \boldsymbol{\phi}(\mathbf{x})
$$
where $\mathbf{E}$ is the macroscopic strain.

In Fourier space, this convolution becomes a simple multiplication:
$$
\hat{\nabla\mathbf{u}}(\mathbf{k}) = - \hat{\boldsymbol{\Gamma}}(\mathbf{k}) : \hat{\boldsymbol{\phi}}(\mathbf{k}) \quad \forall \mathbf{k} \neq \mathbf{0}
$$
The Fourier transform of the Green's operator $\hat{\boldsymbol{\Gamma}}(\mathbf{k})$ is defined as:
$$
\hat{\Gamma}_{ijkl}(\mathbf{k}) = k_j k_l (C^0_{ipkq} k_p k_q)^{-1}
$$
where $k_i$ are the components of the wave vector $\mathbf{k}$.

### Finite Strain Formulation
For finite strains, the formulation is adapted to solve for the velocity gradient $\nabla\mathbf{v}(\mathbf{x})$ in the current configuration. A heterogeneous viscosity tensor $\mathbf{L}^{0}(\mathbf{x})$ is introduced, leading to a similar polarization decomposition:
$$
\boldsymbol{\sigma}(\mathbf{x}) = \mathbf{L}^{0}(\mathbf{x}):\nabla\mathbf{v}(\mathbf{x}) + \boldsymbol{\phi}(\mathbf{x})
$$

The velocity gradient is then updated via:
$$
\nabla\mathbf{v}(\mathbf{x}) = \dot{\mathbf{E}} - FT^{-1}\left(\hat{\boldsymbol{\Gamma}}(\mathbf{k}):\hat{\boldsymbol{\Phi}}(\mathbf{k})\right) \cdot \mathbf{F}^{t,-1}(\mathbf{x})
$$
where $\mathbf{F}$ is the deformation gradient and $\dot{\mathbf{E}}$ is the macroscopic strain rate.

## Material Models

Fierro's EVPFFT solvers support advanced material constitutive models to capture the mechanical response of polycrystalline aggregates.

### 1. Crystal Plasticity (Voce Hardening)
This model captures the anisotropy of single crystals and the evolution of plastic slip on specific crystallographic systems.
*   **Kinematics:** The plastic strain rate $\dot{\boldsymbol{\epsilon}}^{p}$ is sum of shear rates $\dot{\gamma}^{s}$ on all active slip systems $s$:
    $$
    \dot{\boldsymbol{\epsilon}}^{p} = \sum_{s} \mathbf{m}^{s} \dot{\gamma}^{s}
    $$
    where $\mathbf{m}^{s}$ is the Schmid tensor.
*   **Flow Rule:** A power-law relationship relates the shear rate to the resolved shear stress $\tau^s$:
    $$
    \dot{\gamma}^{s} = \dot{\gamma}_{0} \left( \frac{|\tau^s|}{\tau_{c}^{s}} \right)^n \text{sgn}(\tau^s)
    $$
*   **Hardening Law (Voce):** The critical resolved shear stress (CRSS) $\tau_{c}^{s}$ evolves with accumulated shear $\Gamma$:
    $$
    \tau_{c}^{s} = \tau_{0}^{s} + (\tau_{1}^{s} + \theta_{1}^{s}\Gamma) \left( 1 - \exp\left( -\frac{\theta_0^s \Gamma}{\tau_1^s} \right) \right)
    $$
    Parameters:
    *   $\tau_{0}^{s}$: Initial CRSS.
    *   $\tau_{1}^{s}$: Back-extrapolated CRSS.
    *   $\theta_{0}^{s}$: Initial hardening rate.
    *   $\theta_{1}^{s}$: Asymptotic hardening rate.

### 2. J2 Plasticity (von Mises)
Supported in `LS-EVPFFT-J2`, this model assumes isotropic plastic behavior.
*   **Flow Rule:** Based on the J2 invariant of the deviatoric stress tensor.
*   **Hardening:** Isotropic hardening law similar to the Voce form, governed by parameters $\sigma_0, \sigma_1, \theta_0, \theta_1$.

### 3. Elasticity
*   **Anisotropic Elasticity:** Defined by the stiffness tensor $C_{ijkl}$. For cubic crystals, this requires three constants: $C_{11}, C_{12}, C_{44}$.
*   **Isotropic Elasticity:** Defined by Young's Modulus ($E$) and Poisson's Ratio ($\nu$).

## Input File Description

The main input file (e.g., `example_evpfft_standalone_inputfile.txt`) controls the simulation parameters. Below is a breakdown of its sections:

### 1. Grid Dimensions and Phases
```text
1 1 1 12               NPHMX, NMODMX, NTWMMX, NSYSMX
8 8 8                  x-dim, y-dim, z-dim
1                      number of phases (nph)                         
1.  1.  1.             RVE dimensions (delt)
```
*   **NPHMX, etc.:** Maximum number of phases, modes, twinning modes, slip systems (array allocation limits).
*   **x-dim, y-dim, z-dim:** Number of voxels in each direction (must match microstructure file).
*   **RVE dimensions:** Physical length of the RVE sides (e.g., 1.0 1.0 1.0 implies a unit cube).

### 2. Microstructure File
```text
/absolute/path/to/random_microstructure_8x8x8.txt
```
*   Path to the file defining the phase ID and grain orientation (Euler angles) for every voxel.

### 3. Phase Information
```text
*INFORMATION ABOUT PHASE #1
0                          igas(iph)
/absolute/path/to/example_plastic_parameters.txt
/absolute/path/to/example_elastic_parameters.txt
```
*   **igas(iph):** Phase type flag (0 for crystal/solid, 1 for gas/void).
*   **Parameter Files:** Paths to the plastic (slip systems, hardening) and elastic properties files for this phase.

### 4. Test Conditions (Boundary Conditions)
This section defines the macroscopic loading applied to the RVE.
```text
* boundary conditions                                                     
    0       1       1           iudot     |    flag for vel.grad. 
    1       0       1                     |    (0:unknown-1:known)        
    1       1       1                     |            
                                          |                               
   -0.35     0.        0.          udot   |    vel.grad values                  
    0.      -0.35      0.                 |                               
    0.       0.         1.0               |                               
                                          |                               
    1       0        0           iscau    |    flag for Cauchy stress            
            1        0                    |                               
                     0                    |                               
                                          |                               
    0.      0.       0.          scauchy  |    Cauchy stress values              
            0.       0.                   |                               
                     0.                   |                               
```
*   **iudot / udot:** Velocity gradient (strain rate) control. `iudot=1` enforces the corresponding value in `udot`. `iudot=0` means that component is unknown (not enforced).
*   **iscau / scauchy:** Cauchy stress control. `iscau=1` enforces the stress value in `scauchy`.
*   **Mixed Control:** You typically enforce either strain rate (`iudot=1`) OR stress (`iscau=1`) for a given component, but not both.

### 5. Run Conditions and Control
```text
0.00005         eqincr (if ictrl>=0) or tdot (if ictrl=-1)
-1              ictrl (1-6: strain comp, 0: VM eq, -1: tdot)              
*INFORMATION ABOUT RUN CONDITIONS                                         
30              nsteps                                                     
0.0000001       err                                                        
50              itmax                                                      
```
*   **ictrl / eqincr:** Defines the time step or strain increment. `ictrl=-1` means fixed time step `tdot`.
*   **nsteps:** Total number of time steps to run.
*   **err:** Convergence tolerance for the iterative solver.
*   **itmax:** Maximum number of iterations per time step.

### 6. Output Flags
```text
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD update hardening (1) or not (0)?
1               IWTEX write texture (1) or not (0)?
1 30            IWFIELDS,IWSTEP write full fields (1=yes) and interval
```
*   Flags to control restart capability (`IRECOVER`, `ISAVE`), physics updates (`IUPDATE`, `IUPHARD`), and output frequency (`IWTEX`, `IWFIELDS`).
