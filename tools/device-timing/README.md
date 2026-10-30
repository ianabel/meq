# Device timing kit

A self-contained build and timing run of MEQ on a GPU node that has no Claude
Code session on it. It reproduces the `MEASUREMENTS.md` M-242 matrix:
`machine-f-diiid` plus `fixed-f-diiid` at refinement levels 1, 2 and 3, each
run on the host at N threads and on the GPU, in interleaved rounds. It writes
one tarball of logs to bring back. It was written for NERSC Perlmutter
(A100, `sm_80`), and nothing in it is Perlmutter-only except the module names
in `job.slurm` and the defaults below.

**Build it here, run it there.** `tools/device-timing/bundle.sh` (run on the
development machine) packs everything into one tarball:
* MEQ's source as it stands;
* the MFEM MEQ links (`meq-integration`, flattened; also at
  `github.com/ianabel/mfem-meq-integration`, private);
* SUNDIALS 7.5.0;
* netcdf-cxx4 4.3.1, Eigen 3.4.0 and Boost 1.86 (only the parts used);
* the case inputs, including the gmsh mesh `machine-f-diiid` needs, so no
  gmsh is required.

`PROVENANCE.txt` in the bundle says which commit and diff it was made from.
The scripts below are for the GPU machine.

## The four steps

```sh
tar xf meq-device-timing-*.tar.gz && cd meq-device-timing
module load PrgEnv-gnu cudatoolkit/13.2 cray-hdf5 cray-netcdf cmake python
export SUNDIALS_DIR=...    # your CUDA SUNDIALS's install prefix (include/, lib64/);
                          # SUNDIALS_ROOT, as gx's Makefile.perlmutter names it, also works
./preflight.sh            # on a LOGIN node: what is there, what will be built
./build.sh                # on a login node, about 30-60 min; MEQ_JOBS=N to change its 16
./check.sh                # host against GPU, every printed answer: PASS or stop
sbatch job.slurm          # the timing matrix, on one exclusive GPU node
```

### The next time: `./update.sh`, not a new bundle

Once a kit is built, `./update.sh` (on a login node, in place of `build.sh`)
pulls MEQ from `github.com/ianabel/meq` and MFEM from
`github.com/ianabel/mfem-meq-integration`, both `main` by default
(`MEQ_GIT_REF=`/`MFEM_GIT_REF=` for any branch, tag or commit), and rebuilds.
MFEM is rebuilt only when its commit moved, and the other libraries not at
all. It refreshes the scripts here and `cases/` from MEQ's tree, except a
script you edited (`job.slurm`), which it keeps and names. It stops if the
configuration's mesh no longer matches the gmsh mesh in `cases/`, which only a
new bundle (or gmsh) can make. Then `./check.sh` and `sbatch job.slurm` as
before. On the development machine, `mirror.sh` pushes a new MFEM snapshot
to the mirror; MEQ's changes reach the run only once they are pushed.

Then bring back `results-*.tar.gz` (and `preflight.txt`, `check.txt`).
`summarize.py results/` turns it into the tables M-242 prints.

**Send `preflight.txt` back before building if anything in it says MISSING**.
The build fails loudly on a missing piece rather than guessing.

## What the build needs from the machine, and what it supplies itself

| | from the machine | otherwise |
|---|---|---|
| C++17 compiler, CMake ≥ 3.24, Python 3 | `PrgEnv-gnu`, `cmake`, `python` | — |
| CUDA toolkit (`nvcc`) | `cudatoolkit` | — |
| cuDSS 0.8.0 | `CUDSS_DIR` | downloaded from NVIDIA's redistributable, sha256-checked, for the CUDA major `nvcc` reports (login nodes have internet) |
| oneMKL (PARDISO, LAPACK) | `MKLROOT` (an `intel`/`mkl` module) | `python3 -m pip install mkl-devel` into the prefix |
| netCDF-C | `cray-netcdf` (`NETCDF_DIR`): the SERIAL one; gx's `cray-netcdf-hdf5parallel` would pull MPI into a serial code | — |
| netcdf-cxx4 | — | built from the bundle against netCDF-C |
| SUNDIALS with CUDA | `SUNDIALS_DIR`, the install prefix (holding `include/` and `lib/` or `lib64/`), exported: `ENABLE_CUDA=ON` at the same CUDA arch, `SUNDIALS_PRECISION=DOUBLE`, `SUNDIALS_INDEX_SIZE=64`, KINSOL; `build.sh` checks the first two | 7.5.0, built from the bundle |
| Boost.Math, Eigen | `BOOST_ROOT`, `EIGEN3_DIR` | installed from the bundle |

**MKL is what the host arm is measured with.** PARDISO on MKL's GNU OpenMP
layer is MEQ's default host trace solver, and the comparison is only fair
against it. A build that fell back to UMFPACK would be timing a different host
code.

## The timing run

`run.sh` with `job.slurm`'s defaults:

* **cases**: `mf` (`machine-f-diiid`, 4,848 elements, bordered free boundary),
  `ff1`/`ff2`/`ff3` (`fixed-f-diiid` at `RefinementLevels` 1/2/3: 6,245,
  25,487 and 102,840 elements);
* **arms**: `cpu8` (comparable with M-242), `cpu64` (the whole Milan socket),
  `cuda` (one A100);
* **rounds**: 3, interleaved, so drift hits every arm alike; `ff3` takes 2;
* every run writes `--profile`, so the solve leg can be read apart from
  setup and output. The host steps the first Perlmutter run found on one core
  under the Device (`Reconstruct()`, the first assembly, the axis check:
  M-319) are now on the card or threaded (`DEVICE-PLAN.md`, *What is next*);
  the `cuda` arm's profile is where to confirm that held on `ff3`;
* optionally `NSYS=1`: one extra `cuda` run per case under `nsys`, for the
  bus count.

An exclusive node keeps the other three A100s and the 64 cores quiet. The
script records `nvidia-smi`, `lscpu`, the loaded modules and `meq --version`
beside the logs.

## Reading the result

`python3 summarize.py results/` prints, per case and arm, the best round's
wall, setup, solve, output, Newton count, `solve/N` and `psi_ax`, and the
main solve legs. **Compare `solve/N`** where Newton counts differ.
**Distrust a cell whose `psi_ax` disagrees between arms beyond the 7th
digit**: `check.sh` should have stopped the run first.
