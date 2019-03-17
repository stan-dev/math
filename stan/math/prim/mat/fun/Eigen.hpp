#ifndef STAN_MATH_PRIM_MAT_FUN_EIGEN_HPP
#define STAN_MATH_PRIM_MAT_FUN_EIGEN_HPP

#ifdef STAN_OPENCL
#ifdef STAN_OPENCL_CACHE
#include <CL/cl.hpp>
#define EIGEN_MATRIXBASE_PLUGIN "stan/math/prim/mat/matrix_addons.h"
#endif
#endif
#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/src/Core/NumTraits.h>

#endif
