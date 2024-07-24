import platform
nvhpc_ver = '24.5'
cuda_ver = '12.4'
gcc_ver = '13'
llvm_ver = '18'
cmake_ver = '3.27.0'
boost_ver = '1.75.0'
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

# Install and configure AdaptiveCpp:
# Install and configure AdaptiveCpp:
if True:
    Stage0 += boost(version=boost_ver)
    Stage0 += shell(commands=[
        'set -ex',
        # Need this for AdaptiveCpp to find LLVM
        f'ln -sf /usr/lib/llvm-{llvm_ver}/lib/libLLVM-{llvm_ver}.so /usr/lib/llvm-{llvm_ver}/lib/libLLVM.so',
        'git clone --recurse-submodules -b develop https://github.com/AdaptiveCpp/AdaptiveCpp',
        'cd AdaptiveCpp',
        'git submodule update --recursive',
        f'cmake -Bbuild -H.  -DCMAKE_C_COMPILER="$(which clang-{llvm_ver})" -DCMAKE_CXX_COMPILER="$(which clang++-{llvm_ver})" -DCMAKE_INSTALL_PREFIX=/opt/adaptivecpp  -DWITH_CUDA_BACKEND=ON  -DWITH_CPU_BACKEND=ON',
        'cmake --build build --target install -j $(nproc)',
    ])
    Stage0 += environment(variables={
        'PATH':'$PATH:/opt/adaptivecpp/bin',
        'ACPP_APPDB_DIR': '/src/',
    })
