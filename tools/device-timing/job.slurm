#!/bin/bash
# The device timing matrix on one exclusive Perlmutter GPU node. Edit the
# account, then: sbatch job.slurm (from the kit directory).
#
# Exclusive, because the host arms use all 64 cores and the other three A100s
# must be idle for the GPU arm to be the only tenant. About 1 hour with the
# defaults; ff3's arms dominate.
#SBATCH --job-name=meq-device-timing
#SBATCH --account=CHANGE_ME_g
#SBATCH --constraint=gpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --gpus=4
#SBATCH --time=02:00:00
#SBATCH --output=timing-%j.out

module load PrgEnv-gnu cudatoolkit/13.2 cray-hdf5 cray-netcdf cmake python
# The MKL module build.sh used, if any, has to be loaded here too: env.sh
# carries its library path, but a module may also set what MKL reads at run time.

cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
export CUDA_VISIBLE_DEVICES=0

# The gate again on the node the timing runs on: a build that passed on a login
# node's A100 is still checked on the compute node's.
./check.sh || exit 1
./run.sh
