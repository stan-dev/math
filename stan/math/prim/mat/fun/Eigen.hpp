#ifndef STAN_MATH_PRIM_MAT_FUN_EIGEN_HPP
#define STAN_MATH_PRIM_MAT_FUN_EIGEN_HPP

#ifdef EIGEN_MATRIXBASE_PLUGIN
#ifndef EIGEN_STAN_MATRIXBASE_PLUGIN
#error "Stan uses Eigen's EIGEN_MATRIXBASE_PLUGIN macro. To use your own "
"plugin add the eigen_plugin.h file to your plugin."
#endif
#else
#define EIGEN_MATRIXBASE_PLUGIN "stan/math/prim/mat/eigen_plugins.h"
#endif

#ifdef EIGEN_SPARSEMATRIXBASE_PLUGIN
#ifndef EIGEN_STAN_SPARSEMATRIXBASE_PLUGIN
#error "Stan uses Eigen's EIGEN_MATRIXBASE_PLUGIN macro. To use your own "
"plugin add the eigen_plugin.h file to your plugin."
#endif
#else
#define EIGEN_SPARSEMATRIXBASE_PLUGIN "stan/math/prim/mat/eigen_plugins.h"
#endif

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/QR>
#include <Eigen/src/Core/NumTraits.h>

#endif
