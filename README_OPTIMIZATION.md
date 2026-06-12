# JIMWLK solver — optimization & bugfix summary

Base: github.com/hejajama/jimwlk, master @ 6d5d56b.
Modified files: `src/Matrix.{h,cpp}`, `src/Cell.{h,cpp}`, `src/FFT.{h,cpp}`,
`src/Lattice.cpp`, `src/Init.cpp`, `src/main.cpp`, `CMakeLists.txt`.
`optimization.patch` is the full unified diff; the tarball is the complete tree.

## Measured performance (single core; OpenMP adds further scaling on multicore)

| case (simpleLangevin 1, SU(3), MV init)        | original | optimized | speedup |
|------------------------------------------------|----------|-----------|---------|
| 64^2,  Ny=10, 10 steps                          | 0.80 s   | 0.36 s    | 2.2x    |
| 256^2, Ny=20, 6 steps                           | 13.3 s   | 6.4 s     | 2.1x    |
| 700^2, Ny=10, 2 steps                           | 60.8 s   | 33.5 s    | 1.8x    |
| peak RSS at 700^2                               | 2.37 GB  | 1.33 GB   | -44%    |

The old (non-simple) Langevin path speeds up by ~3.0x single-threaded.
All deterministic per-cell loops (the expm-dominated evolution updates, the
adjoint construction, the MV initialization) are OpenMP-parallel; noise
generation stays serial so the RNG sequence per seed is unchanged.

## Validation

Reference = original code with two *identical* test-only patches
(deterministic FFTW plans; text output enabled for MV init). Same seed,
same input, binary full-precision output compared.

- old Langevin path (simpleLangevin 0): **bitwise identical** (SU(2) and SU(3)).
- simple Langevin path, built with `-DLEGACY_SL_SUM=ON`: **bitwise identical**.
- simple Langevin path, default build: differs only by the documented
  factorization $V(\sum_a \xi_a t^a)V^\dagger$ (one re-associated sum):
  max |dU| = 5.6e-16 after 4 steps; dipole amplitude per trajectory agrees to
  ~5e-9 after 50 steps; 6-seed ensemble dipole D(r) identical to all printed
  digits (0.000 sigma pull).
- SU(3) unitarity of evolved Wilson lines: max |UU^dag-1| = 3e-15.
- Note: `-march=native` (now default, disable with `-DENABLE_NATIVE=OFF`)
  enables FMA contraction, which by itself perturbs results at 1 ulp; use
  `-DENABLE_NATIVE=OFF -DDETERMINISTIC_FFT=ON -DLEGACY_SL_SUM=ON` for strict
  bit-compatibility with the original.

## Bugs found and fixed

1. `Matrix::FrobeniusNorm()`: `norm` uninitialized and assigned instead of
   accumulated -> returned |last element|. (Feeds sqrtm/logm convergence,
   used only by the currently disabled Infrared::regulate.)
2. `Matrix::OneNorm()`: column sums overwritten instead of accumulated;
   uninitialized read for 2x2.
3. `Matrix::expm(t,p)`: factor t dropped when ||A||<=1/2 and applied
   per-row in the norm loop (latent; all callers use t=1). Also j1/j2
   uninitialized in logm; VLA replaced by std::vector.
4. `main.cpp` running-coupling kernel S: cos^4(pi x) vs cos^2(pi y) asymmetry
   (must be x<->y symmetric, cf. the fixed-coupling branch) and mass
   regulator applied once instead of squared (S ~ K^2). PHYSICS-CHANGING for
   runningCoupling=1 + simpleLangevin=0 only; does not affect your
   simpleLangevin=1 runs.
5. `Lattice::PrintWilsonLines`: initMethods 1-4 silently wrote no output;
   binary writer hardcoded Nc=3 (wrong for SU(2)); val1 leaked.
6. `Matrix::getElementsText()`: text Wilson lines truncated to 6 significant
   digits (these files are read back by initMethod 10) -> now precision(17).
7. FFT destructor leaked the backward plan; broken+unused fftnComplex removed;
   GSL table leak in logm_pade; delete vs delete[] mismatches in main and Cell.

## Main optimizations

- Matrix: move semantics (the user-declared copy-assignment had suppressed
  them, so every operator temporary was deep-copied), rvalue operator
  overloads, allocation-free conjg(), mult() into preallocated storage,
  in-place axpy addMultiple(), direct-initialized constructors.
- Simple Langevin: $V(\sum_a \xi_a t^a)V^\dagger$ instead of $\sum_a \xi_a(V t^a V^\dagger)$:
  2 triple matrix products + 1 conjugation per site instead of 16 + 8 (the
  original also copied+conjugated U inside the color loop).
- Old Langevin: precomputed sparse adjoint-generator structure (same
  summation order -> bitwise identical), allocation-free computeAdjointU.
- OpenMP on all deterministic per-cell loops (results independent of thread
  count).
- FFT: four quadrant-copy loop nests -> one precomputed shift map.
- Init (MV): rho array allocated once instead of per longitudinal slice;
  expm loop parallelized.
- Cell: matrices by value; the 10 Uy matrices (unequal-rapidity feature)
  allocated lazily -> ~10x fewer startup allocations, main memory saving.
- CMake options: ENABLE_NATIVE (default ON), DETERMINISTIC_FFT,
  LEGACY_SL_SUM (both default OFF).
