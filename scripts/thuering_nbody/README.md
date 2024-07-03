> Comparing a Naive and a Tree-Based N-Body Algorithm using Different Standard SYCL Implementations on Various Hardware

Paper: https://dl.acm.org/doi/10.1145/3624062.3624604
## Building
Follow instructions on the repo (https://github.com/TimThuering/N-Body-Simulation).
As part of the build, you will need to install AdaptiveCpp first.
When building the n-body software, you may find the following `cmake` options useful:
- `-DCMAKE_PREFIX_PATH="/opt/adaptivecpp"` (e.g. AdaptiveCpp location)
- `-DCMAKE_CXX_FLAGS="-I/opt/cuda/targets/x86_64-linux/include"` (e.g. CUDA location)

## Downloading the dataset
Downloading the largest dataset often causes the NASA website (https://ssd.jpl.nasa.gov/tools/sbdb_query.html) to hang.
I have found running `curl` to get the raw json file and then converting to csv to be the best option.

Download `asteroids.json` (this is a ~800 MB file, so may take a while to download):
```bash
curl 'https://ssd-api.jpl.nasa.gov/sbdb_query.api?fields=spkid,full_name,pdes,name,prefix,neo,pha,sats,H,G,M1,M2,K1,K2,PC,diameter,extent,albedo,rot_per,GM,BV,UB,IR,spec_B,spec_T,H_sigma,diameter_sigma,orbit_id,epoch,epoch_mjd,epoch_cal,equinox,e,a,q,i,om,w,ma,ad,n,tp,tp_cal,per,per_y,moid,moid_ld,moid_jup,t_jup,sigma_e,sigma_a,sigma_q,sigma_i,sigma_om,sigma_w,sigma_ma,sigma_ad,sigma_n,sigma_tp,sigma_per,class,producer,data_arc,first_obs,last_obs,n_obs_used,n_del_obs_used,n_dop_obs_used,condition_code,rms,two_body,A1,A1_sigma,A2,A2_sigma,A3,A3_sigma,DT,DT_sigma&full-prec=false&sb-kind=a&www=1' -o asteroids.json
```
Convert to `asteroids.csv`:
```bash
jq -r '
  (.fields | @csv), 
  (.data[] | @csv)
' asteroids.json | tr -d '"' > asteroids.csv
```
Then the dataset can then be preprocessed by the n-body software:
```bash
./dataset_converter/preprocess asteroids.csv -o sim_asteroids.csv
```

## Running
For GPU (default):
```bash
./N_Body_Simulation --file=dataset_converter/sim_asteroids.csv --dt=1h --t_end=1d --vs=1h --vs_dir=./output --algorithm=BarnesHut --theta=0.6
```
To run on CPU add the flag
`--use_gpus=false`.

## Convert to Barnes-Hut format
The following script can be used to convert `sim_asteroids.csv` to a binary format that the main repo n-body software can understand.
```bash
python3 scripts/thuering_nbody/conv_csv.py sim_asteroids.csv sim_asteroids.bin
```
The simulation can then be run (with `D=3`):
```bash
./nbody_d3_gcc -s 25 --workload load scripts/sim_big.bin --print-info --theta 0.6 --precision double
```
This tries to match the above `N_Body_Simulation` as closely as possible:
- In `conv_csv.py` the `.bin` file is saved with `3600` seconds per timestep
- `-s 25` runs for 1 day (including a warm-up timestep)
- `--theta 0.6` the given threshold
- `--precision double` the `N_Body_Simulation` program only uses double precision
