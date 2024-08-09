# n-body mini-apps

![Galaxy collision](./cover_animation.gif)

An implementation of the All-pairs, Barnes-Hut and BVH algorithms for n-body simulations using ISO C++ parallel algorithms.

Two portable execution environments are provided to reproduce:
- Docker
- Mamba

## Running with docker

Pre-requisites: `docker` installed.
Install [HPCCM](https://github.com/NVIDIA/hpc-container-maker):

```shell
$ pip install hpccm
```

Then you can run the samples using:

```shell
# Options
# ./ci/run_docker <toolchain> <algorithm> <workload case> <dim> <precision>
# Example: nvc++ gpu compiler, barnes-hut algorithm, galaxy simulation, 3D, double precision:
$ ./ci/run_docker nvgpu barnes-hut galaxy 3 double
# If you only want to build the binary, you can use:
$ BUILD_ONLY=1 ./ci/run nvgpu barnes-hut galaxy 3 double
# If you only want to run the binary after building it, you can use:
$ RUN_ONLY=1 ./ci/run nvgpu barnes-hut galaxy 3 double
```

Following options available:

- Toolchain: `nvgpu` (`nvc++ -stdpar=gpu`), `nvcpu` (`nvc++ -stdpar=cpu`), `gcc` (Intel TBB), `clang` (Intel TBB), `acpp` ([AdaptiveCpp](https://github.com/AdaptiveCpp/AdaptiveCpp)).
- Algorithm: `all-pairs`, `all-pairs-collapsed`, `barnes-hut`.
- Dimensions: `2` (2D), `3` (3D).
- Precision: `float`, `double`.

When contributing code, you can format your contributions as follows:

```shell
$ ./ci/run_docker fmt
```

## Running with mamba

### Installing

The environment is made portable through mamba/conda.
This must be installed as a prerequisite
e.g. run the Miniforge installer from https://github.com/conda-forge/miniforge.
Then create the `stdpar-bh` environment:
```bash
$ mamba env create -f environment.yaml
```
<!-- `mamba env export --from-history --name stdpar-bh` -->

Other things you might want:
- NVIDIA HPC SDK
<!--- Intel oneAPI Base Toolkit-->


### Building
Use `make` to build the program.
This must be done within the mamba environment:
```bash
$ mamba activate stdpar-bh
```
The number of dimensions can be specified with `D=<dim>` parameter to `make`.
By default `D=2` is used.
These are the available targets:

**CPU**

- `make gcc`
- `make clang`
- `make nvcpp`

**GPU**
- `make gpu` to build for NVIDIA GPUs using `nvc++`

The output will be `./nbody_d<dim>_<target>`.
<!-- makelocalrc -gcc $(which gcc) -gpp $(which g++) -x -d . -->

### Run configuration
When running the `nvcpp` version, it is recommended to use the following environment variables:
```bash
OMP_PLACES=cores OMP_PROC_BIND=close ./nbody_d2_nvcpp -s 5 -n 1000000
```

If you get an error about missing libraries then try running with the following environment variable:
```bash
LD_LIBRARY_PATH=${CONDA_PREFIX}/lib ./nbody_d2_clang -s 5 -n 1000000
```

### Examples
Run Barnes-Hut with $\theta=0$ and compare with all pairs algorithm.
Run 5 steps with 10 bodies.
They should have the same output.
```bash
$ ./nbody_d2_gpu -s 5 -n 10 --print-state --theta 0
$ ./nbody_d2_gpu -s 5 -n 10 --print-state --algorithm all-pairs
```

Run a large Barnes-Hut simulation with 1,000,000 bodies:
```bash
$ ./nbody_d2_gpu -s 5 -n 1000000
```

Generate a similar image to the above GIF:
```bash
$ ./nbody_d2_gpu -s 1000 -n 10000 --save pos --workload galaxy
$ python3 scripts/plotter.py pos --galaxy --gif
```

To find other program arguments:
```bash
$ ./nbody_d2_gpu --help
```
