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

The OpenCL backend can be used for reverse mode AD as well as primitive functions on containers of basic C++ scalar types. Below is the list of functions and distributions that are currently supported.

### Reverse mode

An example of using OpenCL supported Stan Math functions for reverse mode AD is shown below: 

```cpp
using stan::math;

Eigen::Matrix<double, -1, -1> A = Eigen::Matrix<double, -1, -1>::Random(N, M);
Eigen::Matrix<double, -1, -1> B = Eigen::Matrix<double, -1, -1>::Random(M, N);
var_value<matrix_cl<double>> A_cl = to_matrix_cl(A);
var_value<matrix_cl<double>> B_cl = to_matrix_cl(B);
var_value<matrix_cl<double>> C_cl = A_cl * B_cl;
var C_sum = sum(C_cl);
C_sum.grad();
```

#### Supported functions

- acos, acosh
- add, operator+
- add_diag, diag_matrix, diagonal
- asin, asinh
- atan, atanh
- beta, lbeta
- block, col, row, sub_col, sub_row
- cbrt
- ceil, floor, round, trunc
- cholesky_decompose
- cols, dims, num_elements, rows, size
- columns_dot_product, columns_dot_self
- cos, cosh
- crossprod, tcrossprod
- diag_post_multiply, diag_pre_multiply
- digamma
- distance, squared_distance
- dot_product, dot_self
- elt_divide, operator./
- elt_multiply, operator.*
- erf, erfc
- exp, exp2, expm1
- fabs
- hypot
- inv, inv_cloglog, inv_logit, inv_sqrt, inv_square, inv_Phi
- ldexp
- lgamma
- log, log10, log1m, log1p, log2, logit
- log1m_exp, log1p_exp, log1m_inv_logit, log_inv_logit
- log_diff_exp, log_inv_logit_diff
- rows_dot_product, rows_dot_self
- head, tail
- mean
- mdivide_left_tri_low, mdivide_right_tri_low
- multiply, operator*
- multiply_log, lmultiply
- Phi, Phi_approx
- pow
- segment
- sin, sinh
- sqrt
- square
- subtract, operator-
- sum
- tan, tanh
- tgamma
- transpose

#### Supported distributions

- bernoulli_lpmf
- bernoulli_logit_lpmf
- bernoulli_logit_glm_lpmf
- beta_lpdf
- beta_proportion_lpdf
- binomial_lpmf
- categorical_logit_glm_lpmf
- cauchy_lpdf
- chi_square_lpdf
- double_exponential_lpdf
- exp_mod_normal_lpdf
- exponential_lpdf
- frechet_lpdf
- gamma_lpdf
- gumbel_lpdf
- inv_chi_square_lpdf
- inv_gamma_lpdf
- logistic_lpdf
- lognormal_lpdf
- neg_binomial_lpmf
- neg_binomial_2_lpmf
- neg_binomial_2_log_lpmf
- neg_binomial_2_log_glm_lpmf
- normal_lpdf
- normal_id_glm_lpdf
- ordered_logistic_glm_lpmf
- pareto_lpdf
- pareto_type_2_lpdf
- poisson_lpmf
- poisson_log_lpmf
- poisson_log_glm_lpmf
- rayleigh_lpdf
- scaled_inv_chi_square_lpdf
- skew_normal_lpdf
- std_normal_lpdf
- student_t_lpdf
- uniform_lpdf
- weibull_lpdf

### Primitive functions

OpenCL supported primitive functions can be used on `matrix_cl<T>` objects, where T is a primitive built-in C++ type.
An example:

```cpp
using stan::math;

Eigen::Matrix<double, -1, -1> A = Eigen::Matrix<double, -1, -1>::Random(N, M);
Eigen::Matrix<double, -1, -1> B = Eigen::Matrix<double, -1, -1>::Random(M, N);
matrix_cl<double> A_cl = to_matrix_cl(A);
matrix_cl<double> B_cl = to_matrix_cl(B);
matrix_cl<double> C_cl = A_cl * transpose(lgamma(B_cl));
Eigen::Matrix<double, -1, -1> C = from_matrix_cl(C_cl);
```

In addition to the functions listed in the reverse mode list, 
the following functions can be used with `matrix_cl<T>` arguments:

- append_col, append_row
- gp_exp_quad_cov
- trace
- trigamma
