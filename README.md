# Barnes-Hut mini-app

![Galaxy collision](./cover_animation.gif)

An implementation of the Barnes-Hut algorithm using C++17 parallel algorithms.

## Installing
<!--
Download the submodules:
```bash
$ git submodule init
$ git submodule update
```
-->
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


## Building
Use `make` to build the program.
This must be done within the mamba environment:
```bash
$ mamba activate stdpar-bh
```
These are the available targets:

**CPU**

- `make gcc`
- `make clang`
- `make nvcpp`

**GPU**
- `make gpu` to build for NVIDIA GPUs using `nvc++`

The output will be `./nbody_<target>`.
<!-- makelocalrc -gcc $(which gcc) -gpp $(which g++) -x -d . -->

## Run configuration
When running the `nvcpp` version, it is recommended to use the following environment variables:
```bash
OMP_PLACES=cores OMP_PROC_BIND=close ./nbody_nvcpp -s 5 -n 1000000
```

If you get an error about missing libraries then try running with the following environment variable:
```bash
LD_LIBRARY_PATH=${CONDA_PREFIX}/lib ./nbody_clang -s 5 -n 1000000
```

## Examples
Run Barnes-Hut with $\theta=0$ and compare with all pairs algorithm.
Run 5 steps with 10 bodies.
They should have the same output.
```bash
$ ./nbody_gpu -s 5 -n 10 --print-state --theta 0
$ ./nbody_gpu -s 5 -n 10 --print-state --all-pairs
```

Run a large Barnes-Hut simulation with 1,000,000 bodies:
```bash
$ ./nbody_gpu -s 5 -n 1000000
```

Generate a similar image to the above GIF:
```bash
$ ./nbody_gpu -s 1000 -n 10000 --save-pos --galaxy
$ python3 plotter.py pos --galaxy --gif
```

To find other program arguments:
```bash
$ ./nbody_gpu --help
```
