#ifndef STAN_MATH_PRIM_FUN_EIGEN_HPP
#define STAN_MATH_PRIM_FUN_EIGEN_HPP

#ifdef EIGEN_DENSEBASE_PLUGIN
#ifndef EIGEN_STAN_DENSEBASE_PLUGIN
#error "Stan uses Eigen's EIGENDENSEBASE_PLUGIN macro. To use your own "
"plugin add the eigen_plugin.h file to your plugin."
#endif
#else
#define EIGEN_DENSEBASE_PLUGIN "stan/math/prim/eigen_plugins.h"
#endif

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/QR>
#include <Eigen/src/Core/NumTraits.h>
#include <Eigen/SVD>

    namespace Eigen {

  /**
   * Traits specialization for Eigen binary operations for `int`
   * and `double` arguments.
   *
   * @tparam BinaryOp type of binary operation for which traits are
   * defined
   */
  template <typename BinaryOp>
  struct ScalarBinaryOpTraits<int, double, BinaryOp> {
    using ReturnType = double;
  };

  /**
   * Traits specialization for Eigen binary operations for `double`
   * and `int` arguments.
   *
   * @tparam BinaryOp type of binary operation for which traits are
   * defined
   */
  template <typename BinaryOp>
  struct ScalarBinaryOpTraits<double, int, BinaryOp> {
    using ReturnType = double;
  };

  /**
   * Traits specialization for Eigen binary operations for `int`
   * and complex `double` arguments.
   *
   * @tparam BinaryOp type of binary operation for which traits are
   * defined
   */
  template <typename BinaryOp>
  struct ScalarBinaryOpTraits<int, std::complex<double>, BinaryOp> {
    using ReturnType = std::complex<double>;
  };

  /**
   * Traits specialization for Eigen binary operations for complex
   * `double` and `int` arguments.
   *
   * @tparam BinaryOp type of binary operation for which traits are
   * defined
   */
  template <typename BinaryOp>
  struct ScalarBinaryOpTraits<std::complex<double>, int, BinaryOp> {
    using ReturnType = std::complex<double>;
  };

}  // namespace Eigen

#endif
