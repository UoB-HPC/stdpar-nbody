import platform
nvhpc_ver = '24.5'
cuda_ver = '12.4'
rocm_ver = '6.1.3'
gcc_ver = '13'
llvm_ver = '18'
gpu = USERARG.get('gpu', 'nv')
if gpu == 'amd' and rocm_ver == '6.1.3':
    # Need to work around lack of support for gcc 13 in rocm-stdpar
    # and different LLVM ABI in LLVM 18 than 18 for AdaptiveCpp.
    gcc_ver = '12'
    llvm_ver = '17'
cmake_ver = '3.27.0'
ubuntu_ver = '22.04'
arch = hpccm.config.get_cpu_architecture()

# NVIDIA HPC SDK base container image
img = f'nvcr.io/nvidia/nvhpc:{nvhpc_ver}-devel-cuda{cuda_ver}-ubuntu{ubuntu_ver}'
Stage0 += baseimage(image = img)

# Install compilers: GCC, LLVM/Clang...
Stage0 += gnu(version=gcc_ver, extra_repository=True)
Stage0 += llvm(version=llvm_ver, upstream=True, extra_tools=True, toolset=True, _trunk_version='19')

# Install other relevant packages: Intel TBB for libstdc++ parallel algorithms, libnuma for NUMA controls,
# python for visualization, etc.
Stage0 += cmake(eula=True, version=cmake_ver)
Stage0 += packages(ospackages=[
    'libtbb-dev', 'libnuma-dev', 'numactl', 'git', 'make', 'bc', 'curl', 'nginx', 'build-essential', 'wget',
    'python3', 'python3-pip', 'python-is-python3', 'python3-setuptools', 'python3-dev',
    'libboost-all-dev', ' libfmt-dev',
])

# Setup miscelanous things:
Stage0 += shell(commands=[
    'set -ex',  # Exit on first error and debug output

    # Workaround docker runtime bug that fails to link libnvidia-ml as .so:
    f'ln -s /usr/lib/{arch}-linux-gnu/libnvidia-ml.so.1 /usr/lib/{arch}-linux-gnu/libnvidia-ml.so',
    
    # C++23 has a feature (https://wg21.link/p2408) that enables iterators from random access ranges to be used with the parallel algorithms.
    # This feature can be backported to older C++ versions.
    # GCC does not implement this feature yet, so when using a `std::views::iota().begin()` iterator with parallel algorithms, libstdc++ runs the algorithm sequentially.
    # Clang and NVC++ already implement this and don't have this issue.
    # GCC bug tracking this is: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=110512
    # The following command patches libstdc++ to workaround the lack of this feature in GCC:
    f"sed -i 's@std::is_same<typename std::iterator_traits<_IteratorType>::iterator_category, std::random_access_iterator_tag>@std::integral_constant<bool, std::random_access_iterator<_IteratorType>>@g' /usr/include/c++/{gcc_ver}/pstl/execution_impl.h",
    
    # Update the HPC SDK toolchain configuration, such that it automatically picks the latest installed GCC:
    f'cd /opt/nvidia/hpc_sdk/Linux_{arch}/{nvhpc_ver}/compilers/bin/',
    'makelocalrc -d . -x .',

    # Install required python packages for the notebooks
    'pip install --upgrade pip',
    'pip install numpy matplotlib gdown jupyterlab ipywidgets pandas seaborn conan jupyterlab-nvidia-nsight',
])

# Install Intel OpenCL and ROCm on x86_64 builds:
acpp_flags=''
if gpu == 'amd' and arch == 'x86_64':
    # Install Intel OpenCL:
    Stage0 += shell(commands=[
        'set -ex',
        'mkdir -p /var/tmp',
        'cd /var/tmp',
        'wget https://github.com/intel/intel-graphics-compiler/releases/download/igc-1.0.16695.4/intel-igc-core_1.0.16695.4_amd64.deb',
        'wget https://github.com/intel/intel-graphics-compiler/releases/download/igc-1.0.16695.4/intel-igc-opencl_1.0.16695.4_amd64.deb',
        'wget https://github.com/intel/compute-runtime/releases/download/24.17.29377.6/intel-level-zero-gpu-dbgsym_1.3.29377.6_amd64.ddeb',
        'wget https://github.com/intel/compute-runtime/releases/download/24.17.29377.6/intel-level-zero-gpu_1.3.29377.6_amd64.deb',
        'wget https://github.com/intel/compute-runtime/releases/download/24.17.29377.6/intel-opencl-icd-dbgsym_24.17.29377.6_amd64.ddeb',
        'wget https://github.com/intel/compute-runtime/releases/download/24.17.29377.6/intel-opencl-icd_24.17.29377.6_amd64.deb',
        'wget https://github.com/intel/compute-runtime/releases/download/24.17.29377.6/libigdgmm12_22.3.19_amd64.deb',
        'wget https://github.com/oneapi-src/level-zero/releases/download/v1.13.5/level-zero-devel_1.13.5+u22.04_amd64.deb',
        'wget https://github.com/oneapi-src/level-zero/releases/download/v1.13.5/level-zero_1.13.5+u22.04_amd64.deb',
        'dpkg -i *.deb',
        'cd -',
        'rm -rf /var/tmp',
    ])
    # Install ROCm
    Stage0 += shell(commands=[
        'set -ex',
        'wget -q -O - https://repo.radeon.com/rocm/rocm.gpg.key | apt-key add -',
        f'echo "deb [arch=amd64] https://repo.radeon.com/rocm/apt/{rocm_ver} focal main" | tee /etc/apt/sources.list.d/rocm.list',
        'printf "Package: *\\nPin: release o=repo.radeon.com\\nPin-Priority: 600" | tee /etc/apt/preferences.d/rocm-pin-600',
        'apt-get update',
        'apt-get install -y rocm-dev rocthrust-dev',
    ])
    #acpp_flags += '-DWITH_OPENCL_BACKEND=ON'
    Stage0 += environment(variables={'AMDCLANG': f'/opt/rocm-{rocm_ver}/lib/llvm/bin/clang++'})
    acpp_flags += '-DWITH_ROCM_BACKEND=ON'

# Install and configure AdaptiveCpp:
if True:
    Stage0 += shell(commands=[
        'set -ex',
        'git clone -b nbody https://github.com/AdaptiveCpp/AdaptiveCpp',
        'cd AdaptiveCpp',
        'git submodule update --recursive',
        f'cmake -Bbuild -H.  -DCMAKE_C_COMPILER="$(which clang-{llvm_ver})" -DCMAKE_CXX_COMPILER="$(which clang++-{llvm_ver})" -DCMAKE_INSTALL_PREFIX=/opt/adaptivecpp -DWITH_CUDA_BACKEND=ON {acpp_flags}',
        'cmake --build build --target install -j $(nproc)',
    ])
    Stage0 += environment(variables={
        'PATH':'$PATH:/opt/adaptivecpp/bin',
        'ACPP_APPDB_DIR': '/src/',
    })

Stage0 += environment(variables={'MPLCONFIGDIR':'/src/.mpl'})

if False:
    # Get for performance comparsions: https://github.com/TimThuering/N-Body-Simulation
    # Comparing a Naive and a Tree-Based N-Body Algorithm using Different Standard SYCL Implementations on Various Hardware
    # Paper: https://dl.acm.org/doi/10.1145/3624062.3624604
    Stage0 += shell(commands=[
        'set -ex',
        'cd /usr/local',
        'git clone https://github.com/TimThuering/N-Body-Simulation.git',
        'cd -',
    ])

    # Largest NASA dataset from https://ssd.jpl.nasa.gov/tools/sbdb_query.html
    Stage0 += shell(commands=[
        'set -ex',
        'mkdir -p datasets',
        # Download dataset
        "curl -s 'https://ssd-api.jpl.nasa.gov/sbdb_query.api?fields=spkid,full_name,pdes,name,prefix,neo,pha,sats,H,G,M1,M2,K1,K2,PC,diameter,extent,albedo,rot_per,GM,BV,UB,IR,spec_B,spec_T,H_sigma,diameter_sigma,orbit_id,epoch,epoch_mjd,epoch_cal,equinox,e,a,q,i,om,w,ma,ad,n,tp,tp_cal,per,per_y,moid,moid_ld,moid_jup,t_jup,sigma_e,sigma_a,sigma_q,sigma_i,sigma_om,sigma_w,sigma_ma,sigma_ad,sigma_n,sigma_tp,sigma_per,class,producer,data_arc,first_obs,last_obs,n_obs_used,n_del_obs_used,n_dop_obs_used,condition_code,rms,two_body,A1,A1_sigma,A2,A2_sigma,A3,A3_sigma,DT,DT_sigma&full-prec=false&sb-kind=a&www=1' -o /datasets/asteroids.json",
    ])
    Stage0 += shell(commands=[
        "jq -r '\\n  (.fields | @csv),\\n  (.data[] | @csv)\\n' /datasets/asteroids.json | tr -d '\"' > /datasets/asteroids.csv",
        '/usr/local/N-Body-Simulation/dataset_converter/preprocess /datasets/asteroids.csv -o /datasets/sim_asteroids.csv',
    ])
