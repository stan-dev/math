# OpenCL CPU/GPU Support {#opencl_support}

[OpenCL](https://www.khronos.org/opencl/) is an open-source framework for writing programs that utilize a platform with heterogeneous hardware. Stan uses OpenCL to design the GPU routines for the Cholesky Decomposition and it's derivative. Other routines will be available in the future. These routines are suitable for programs which require solving large `NxM` matrices (`N>600`) such as algorithms that utilize large covariance matrices.

# Requirements

Users must have suitable hardware (e.g. Nvidia or AMD gpu), valid OpenCL driver, SDK and a suitable C/C++ compiler installed on their computer.

# Installation

## Linux
The following guide is for Ubuntu, but it should be similar for any other Linux distribution. You should have the GNU compiler suite or clang compiler installed beforehand.

Install the Nvidia CUDA toolkit and clinfo tool if you have an Nvidia GPU
```bash
apt update
apt install nvidia-cuda-toolkit clinfo
```

Those with AMD devices can install the OpenCL driver available through

```bash
apt install -y libclc-amdgcn mesa-opencl-icd clinfo
```

If your device is not supported by the current drivers available you can try Paulo Miguel [PPA](https://laanwj.github.io/2016/05/06/opencl-ubuntu1604.html)

```bash
add-apt-repository ppa:paulo-miguel-dias/mesa
apt-get update
apt-get install libclc-amdgcn mesa-opencl-icd
```

## MacOS
Mac's should already have the OpenCL driver installed if you have the appropriate hardware.

Note that if you are building on a mac laptop you may not have a GPU device. You can still use the OpenCL routines for parallelization on your CPU.

## Windows
Install the latest [Rtools](https://cran.r-project.org/bin/windows/Rtools/) suite if you don't already have it. During the installation make sure that the 64 bit toolchain is installed. You also need to verify that you have the System Enviroment variable `Path` updated to include the path to the g++ compiler (`<Rtools installation path>\mingw_64\bin`).

If you have a Nvidia card, install the latest [Nvidia CUDA toolkit](https://developer.nvidia.com/cuda-toolkit).
AMD users should use [AMD APP SDK](http://amd-dev.wpengine.netdna-cdn.com/app-sdk/installers/APPSDKInstaller/3.0.130.135-GA/full/AMD-APP-SDKInstaller-v3.0.130.135-GA-windows-F-x64.exe).

Users can check that their installation is valid by downloading and running clinfo.  
- [download clinfo.exe](https://ci.appveyor.com/api/projects/oblomov/clinfo/artifacts/clinfo.exe?job=platform%3a+x64)

## Setting up the Math Library to run on a GPU

To turn on GPU computation:

1. Check and record what device and platform you would like to use with [clinfo](https://github.com/Oblomov/clinfo); you will the platform and device index such as the printout below

```bash
clinfo -l
# Platform #0: Clover
# Platform #1: Portable Computing Language
#  `-- Device #0: pthread-AMD Ryzen Threadripper 2950X 16-Core Processor
# Platform #2: NVIDIA CUDA
#  +-- Device #0: TITAN Xp
#  `-- Device #1: GeForce GTX 1080 Ti
```

2. In the top level of the math library, open a text file called make/local. if it does not exist, create one.
3. Add these lines to the make/local file:

```bash
STAN_OPENCL=true
OPENCL_DEVICE_ID=${CHOSEN_INDEX}
OPENCL_PLATFORM_ID=${CHOSEN_INDEX}
```

where the user will replace ${CHOSEN_INDEX} with the index of the device and platform they would like to use. In most cases these two will be 0. If you are using Windows append the following lines at the end of the make/local file in order to link with the appropriate OpenCL library:

* Nvidia
```make
CC = g++
LDFLAGS_OPENCL= -L"$(CUDA_PATH)\lib\x64" -lOpenCL
```

* AMD
```make
CC = g++
LDFLAGS_OPENCL= -L"$(AMDAPPSDKROOT)lib\x86_64" -lOpenCL
```

## Running Tests with OpenCL

Once you have done the above step, `runTests.py` should execute with the GPU enabled. All tests will match the phrase `*_opencl_*` and tests can be filtered such as

```bash
./runTests.py test/unit -f opencl
```

## Using the OpenCL backend

We currently have support for the following methods

- bernoulli_logit_glm
- cholesky_decompose
- categorical_logit_glm
- gp_exp_quad_cov
- mdivide_right_tri
- mdivide_left_tri
- multiplication
- neg_binomial_2_log_glm
- normal_id_glm
- ordered_logistic_glm
- poisson_log_glm

In cmdstan, an example model is provided in `examples/GP/`, which uses OpenCL Cholesky decomposition. You can check if your OpenCL configuration works by trying to build it.
* Linux and MacOS: `make examples/GP/gp3`
* Windows: `make examples/GP/gp3.exe`
