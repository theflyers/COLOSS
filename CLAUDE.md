# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**COLOSS** (Complex-scaled Optical and couLOmb Scattering Solver) is a Fortran/C++ scientific computing application for nuclear physics scattering calculations. It solves nuclear scattering problems using complex scaling techniques to compute S-matrices and cross-sections for nuclear reactions (e.g., neutron-nucleus, alpha-nucleus scattering).

**Key Technologies:**
- Fortran (primary computational code)
- C++ (high-performance Coulomb wave function library in `adyo_v1_0/`)
- LAPACK for linear algebra
- OpenMP for parallelization

## Build System

### Prerequisites
- GCC/gfortran compiler
- LAPACK library
- Standard C++11 compiler

### Compilation Commands

**Important:** Edit the Makefile first to set your LAPACK path:
```makefile
LIB = -L/path/to/your/lapack -llapack
```

**Build from scratch:**
```bash
make clean         # Remove all build artifacts
make              # Compile everything (builds adyo_v1_0/, then main program)
```

The build process:
1. Compiles C++ library `adyo_v1_0/libcwf_cpp.a` (Coulomb wave functions)
2. Compiles all Fortran modules into `liballmodules.a`
3. Links final `COLOSS` executable

**Build flags:**
- `-fopenmp`: OpenMP parallelization enabled
- `-g`: Debug symbols included
- `-cpp`: C preprocessor for Fortran

## Running COLOSS

**Basic execution:**
```bash
./COLOSS < input_file
```

**Test cases available in `test/` directory:**
```bash
cd test/
../COLOSS < n40Ca.in      # Neutron on 40Ca at 20 MeV
../COLOSS < alpha40Ca.in  # Alpha on 40Ca
../COLOSS < d93Nb.in      # Deuteron on 93Nb
../COLOSS < alpha28Si.in  # Alpha on 28Si
../COLOSS < 6Li208Pb.in   # 6Li on 208Pb
../COLOSS < readtest.in   # External potential from pot.dat
../COLOSS < yama.in       # Yamaguchi potential test
```

**Performance timing:**
```bash
./timer.sh  # Runs 100 iterations to measure performance
```

## Code Architecture

### Main Program Flow (COLOSS.F)

The execution follows this pipeline:
1. `read_input()` - Parse 4 namelists from stdin
2. `generate_channels()` - Create L-S-J quantum number combinations
3. `init_laguerre_mesh()` - Generate complex-scaled Lagrange-Laguerre basis mesh
4. `get_pot_para()` - Initialize optical potential parameters
5. `basis_func()` - Compute basis functions
6. `init_coul()` - Initialize Coulomb wave functions (calls C++ library)
7. **Solve scattering** (method-dependent):
   - Method 1: `solve_scatt()` - Linear equation approach
   - Method 2: `solve_bound()` + `solve_scatt_green()` - Green's function approach
8. `xsec()` - Calculate differential cross-sections
9. `outinfo()` - Write results to fort.X files

### Core Module Organization

**Foundation modules** (define types and constants):
- `precision.F90` - SP/DP types, machine epsilon, pi, imaginary unit
- `constants.F90` - Physical constants (ℏc, α, amu, e²)
- `system.f` - Reaction system parameters (masses, charges, energy, J ranges)
- `mesh.f` - Mesh arrays and complex scaling angle
- `pot_class.f` - Optical model potential container

**Computational modules:**
- `generate_laguerre.f` - Lagrange-Laguerre basis functions, interfaces with C++ Coulomb library
- `matrix_element.f` - Construct Hamiltonian matrix elements (kinetic + potential)
- `solve_eigen.f` - LAPACK eigenvalue solver wrapper (zggev)
- `bound.f` - Bound state calculation for Green's function method
- `scatt.f` - Scattering amplitude and S-matrix computation
- `rot_potential.f` - Complex rotation operations

**Auxiliary modules:**
- `channels.f` - Generate L-S-J channel indices
- `input.f` - Namelist I/O handling
- `readpot.f` - Import external potential from `pot.dat`
- `clebsch.f` - Clebsch-Gordan coefficients
- `spharm.f` - Spherical harmonics and Legendre polynomials
- `gauss.F`, `gauss_mesh.f90` - Gauss-Legendre quadrature
- `coulcc.f`, `coul90.f` - Fortran Coulomb wave functions (alternatives)
- `yamaguchi.f` - Test potential for validation

### Key Algorithmic Details

**Complex Scaling:**
- Rotates coordinates by angle `ctheta` to transform oscillatory boundary conditions into exponentially decaying ones
- Two rotation modes controlled by `backrot`:
  - `.false.`: Rotate potential directly
  - `.true.`: Backward rotate basis functions

**Basis Functions:**
- Generalized Laguerre polynomials L_n^α(x)
- Lagrange function interpolation for mesh points
- Parameter `alpha` controls polynomial shape
- Typical convergence: 50-100 basis states (`nr`)

**Numerical Integration:**
- Default: Approximated Lagrange function integration
- Optional: Gauss-Legendre quadrature (`matgauss=t`, `bgauss=t`)

**Coulomb Wave Functions:**
- Two implementations: `cwftype=1` (COULCC) or `cwftype=2` (cwfcomplex from C++)
- C++ library `adyo_v1_0/` provides high-performance complex Coulomb functions

## Input Format (Namelists)

COLOSS uses Fortran namelists. Four required namelists:

### &general
Key parameters:
- `nr` - Number of basis functions (20-100+)
- `alpha` - Laguerre polynomial parameter
- `Rmax` - Maximum mesh radius (fm)
- `ctheta` - Complex scaling angle (degrees, typically 5-10)
- `method` - Solution method (1=linear equation, 2=Green's function)
- `backrot` - Rotation target (.true.=basis, .false.=potential)
- `thetah`, `thetamax` - Angular distribution output parameters

### &system
Defines the nuclear reaction:
- `zp`, `zt` - Projectile/target charge
- `massp`, `masst` - Projectile/target mass number
- `namep`, `namet` - Names (for output)
- `jmin`, `jmax` - Total angular momentum range
- `sp` - Projectile spin
- `elab` - Lab frame kinetic energy (MeV)

### &pot
Optical Model Potential (Woods-Saxon form):
- Volume terms: `vv, rv, av` (real); `wv, rw, aw` (imaginary)
- Surface terms: `vs, rvs, avs` (real); `ws, rws, aws` (imaginary)
- Spin-orbit: `vsov, rsov, asov` (real); `vsow, rsow, asow` (imaginary)
- Coulomb radius: `rc`

### &nonlocalpot
- `nonlocal` - Enable Perey-Buck nonlocal form
- `nlbeta` - Nonlocality parameter

## Output Files

Results written to Fortran unit numbers:
- `fort.1` - Copy of input file
- `fort.2` - Channel index table (L, S, J)
- `fort.10` - Mesh points and weights
- `fort.60` - S-matrices for each channel
- `fort.61` - Nuclear scattering amplitudes
- `fort.67` - Differential cross-section vs angle

## Development Workflow

### Modifying Computational Code

When changing core algorithms:
1. **Module dependencies:** Check module `use` statements - `precision` and `constants` are used everywhere
2. **Linear algebra:** LAPACK routines are in `solve_eigen.f` (zggev, zgesv, zgeev)
3. **Matrix elements:** Modify `matrix_element.f` for Hamiltonian changes
4. **Scattering calculation:** Edit `scatt.f` for S-matrix/amplitude changes

### Adding New Potentials

Two approaches:
1. **Analytical potential:** Add module similar to `yamaguchi.f`, integrate in main loop
2. **Numerical potential:** Use `readpot.f` interface, prepare `pot.dat` file

Set `readinpot=t` and `backrot=t` to read external potential.

### Debugging Tips

**Check convergence:**
- Increase `nr` (basis size) and verify S-matrix convergence
- Adjust `Rmax` to ensure wave function decays at boundary
- Vary `ctheta` within 5-15 degrees

**Verify physical results:**
- S-matrix should satisfy unitarity: |S|² ≤ 1
- Elastic cross-section should match experimental data
- Check reaction cross-section positivity

**Common issues:**
- LAPACK errors: Check Hamiltonian matrix conditioning
- Coulomb function failures: Verify Sommerfeld parameter range
- Non-converged results: Increase basis size or adjust mesh parameters

## Testing

The `test/` directory contains validated benchmark cases:
- **n40Ca.in**: Reference case for neutron scattering, quick runtime
- **alpha40Ca.in**: Charged particle scattering with Coulomb
- **d93Nb.in**: Used in performance timing
- **readtest.in**: Tests external potential reading

Run test cases after modifications to verify no regressions.

## Important Implementation Notes

### Fortran-C++ Interface
- `adyo_v1_0/cwf_cpp.f90` bridges Fortran to C++ Coulomb functions
- Uses `iso_c_binding` for interoperability
- C++ library must be compiled first (Makefile handles this)

### OpenMP Parallelization
- Enabled via `-fopenmp` flag
- Parallel regions in matrix construction loops
- Performance scales with number of channels and basis size

### Version Information
- Git commit hash and date embedded in binary via preprocessor macros
- See Makefile lines 10-17: `VERDATE`, `VERREV`, `COMPDATE`

## References

Published paper: https://www.sciencedirect.com/science/article/abs/pii/S0010465525000712

Input file generator UI: `coloss-input-generator.html` (check repository root or documentation)
