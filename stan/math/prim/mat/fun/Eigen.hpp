#ifndef STAN_MATH_PRIM_MAT_FUN_EIGEN_HPP
#define STAN_MATH_PRIM_MAT_FUN_EIGEN_HPP

#ifdef STAN_OPENCL
#include <CL/cl.hpp>
#define EIGEN_MATRIXBASE_PLUGIN "stan/math/prim/mat/matrix_addons.h"
#endif
#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/src/Core/NumTraits.h>

#endif
