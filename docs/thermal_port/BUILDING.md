# Workstation build and baseline validation

## Hardware and dependencies

Tested on Fedora 44, Ryzen 9 9950X3D2 (16 cores / 32 threads), with two NVIDIA
RTX PRO 6000 Blackwell Max-Q GPUs. Driver: 610.57.04. Single-GPU tests use physical
GPU 1 while GPU 0 serves the desktop and another application.

Load the available modules:

```sh
source /usr/share/lmod/lmod/init/bash
module load gcc/14.4.0 cuda/12.9.2 openmpi/5.0.10 cmake/3.31.12 bison/3.8.2
export NVCC_WRAPPER_DEFAULT_COMPILER="$(command -v g++)"
```

The automation environment inherits malformed exported Bash functions. The local
`work/run-env.py` wrapper removes only inherited `BASH_FUNC_*` entries and reloads
these modules. It does not change system or user shell configuration.

| Dependency | Pinned source |
|---|---|
| Kokkos 4.7.04 | `82799e4577568f9666bde36265ac15d78da3e6c8` |
| MATAR (current branch baseline) | `9926149c7819e112420e51e30df3b62ea85cbe09` + `PatchMatarKokkos47.cmake` |
| MATAR (student baseline) | `8142fd6e47615be0d6b0401dbc4118eb8ce87dd8` |
| Scotch 7.0.4 | `82ec87f558f4acb7ccb69a079f531be380504c92` + existing ELEMENTS Scotch CMake patch |

Kokkos 4.7 is the first series supporting `BLACKWELL120`; it retains C++17.
See the [official architecture guide](https://kokkos.org/kokkos-core-wiki/get-started/configuration-guide.html).
The CPU compatibility patch ensures MATAR uses the same allocation for host
and device views when the backend is host-accessible. It preserves separate
host/device allocations and synchronization for CUDA.

## CPU

From the Fierro repository root:

```sh
cmake -S apps/multiphysics -B build-cpu -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
  -DFIERRO_ENABLE_OPENMP=ON -DFIERRO_BUILD_TESTS=ON
cmake --build build-cpu -j 8
OMP_NUM_THREADS=4 OMP_PROC_BIND=spread OMP_PLACES=cores \
  ctest --test-dir build-cpu --output-on-failure
```

## CUDA

Make Kokkos 4.7.04's `bin/nvcc_wrapper` available before configuring the compiler
(a source archive is sufficient). With its absolute location in `KOKKOS_SOURCE`:

```sh
cmake -S apps/multiphysics -B build-cuda -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_CXX_COMPILER="$KOKKOS_SOURCE/bin/nvcc_wrapper" \
  -DFIERRO_ENABLE_CUDA=ON -DKokkos_ARCH_BLACKWELL120=ON \
  -DFIERRO_BUILD_TESTS=ON
cmake --build build-cuda -j 8
CUDA_VISIBLE_DEVICES=1 ctest --test-dir build-cuda --output-on-failure
```

All builds use MPI. Test each binary directly and with `mpirun --bind-to none
-np 2`. Two ranks sharing GPU 1 are valid MPI/GPU tests, but are distinct from
two-GPU scaling measurements. Sandbox restrictions hide GPU device nodes and
deny MPI sockets; approved execution outside the sandbox is required here.

Normal FetchContent population was configured in `work/builds/main-default-cmake`:
the downloaded Kokkos and MATAR SHAs match the table and all three MATAR
constructors were patched automatically. Local baseline builds use explicit
`FETCHCONTENT_SOURCE_DIR_*` paths to avoid repeated downloads. Such overrides
bypass FetchContent's patch step, so the baseline harness uses separately patched
MATAR source copies; the original submodules remain unchanged.

## Recorded baseline results

| Repository / test | CPU direct | CPU 2 ranks | CUDA direct | CUDA 2 ranks |
|---|---|---|---|---|
| Current main: SGTM cooling cube | Pass | Pass | Pass | Pass |
| Current main: Sod_X | Pass | Pre-existing reference discrepancy | Pass | Same pre-existing reference discrepancy |
| Student: Sod_X (baseline repairs) | Pass | Pass | Pass | Pass |
| Student: SGTM cooling cube | Missing state output | Missing state output | Host read of GPU memory aborts | Host read of GPU memory aborts |

Current-main cooling comparisons cover 8000 material-state values with maximum
absolute difference about `1.02e-10`. This is regression evidence, not yet analytic
heat-transfer verification. Student thermal failures are not accepted as reference
solutions for the port.

Complete commands, solver logs, outputs and JSON results are under `work/runs`.
`work/regression.py` checks return status, output existence/schema, finite values,
global element keys and reference tolerances. It handles both state-file formats,
preserves failed runs, limits log growth, and terminates timed-out process groups.
