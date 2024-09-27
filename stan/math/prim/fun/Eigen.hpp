#ifndef STAN_MATH_PRIM_FUN_EIGEN_HPP
#define STAN_MATH_PRIM_FUN_EIGEN_HPP

#ifdef EIGEN_MATRIXBASE_PLUGIN
#ifndef EIGEN_STAN_MATRIXBASE_PLUGIN
#error "Stan uses Eigen's EIGEN_MATRIXBASE_PLUGIN macro. To use your own "
"plugin add the eigen_plugin.h file to your plugin."
#endif
#else
#define EIGEN_MATRIXBASE_PLUGIN "stan/math/prim/eigen_plugins.h"
#endif

#ifdef EIGEN_ARRAYBASE_PLUGIN
#ifndef EIGEN_STAN_ARRAYBASE_PLUGIN
#error "Stan uses Eigen's EIGEN_ARRAYBASE_PLUGIN macro. To use your own "
    "plugin add the eigen_plugin.h file to your plugin."
#endif
#else
#define EIGEN_ARRAYBASE_PLUGIN "stan/math/prim/eigen_plugins.h"
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

//   /**
//    * Traits specialization for Eigen binary operations for `int`
//    * and complex `double` arguments.
//    *
//    * @tparam BinaryOp type of binary operation for which traits are
//    * defined
//    */
//   template <typename BinaryOp>
//   struct ScalarBinaryOpTraits<int, stan::math::complex<double>, BinaryOp> {
//     using ReturnType = stan::math::complex<double>;
//   };

//   /**
//    * Traits specialization for Eigen binary operations for complex
//    * `double` and `int` arguments.
//    *
//    * @tparam BinaryOp type of binary operation for which traits are
//    * defined
//    */
//   template <typename BinaryOp>
//   struct ScalarBinaryOpTraits<stan::math::complex<double>, int, BinaryOp> {
//     using ReturnType = stan::math::complex<double>;
//   };

// template <typename V>
// struct NumTraits<stan::math::complex<V>> {
//   using Real = V;
//   using NonInteger = stan::math::complex<V>;
//   using Literal = stan::math::complex<V>;
//   using Nested = stan::math::complex<V>;
//   enum {
//     IsComplex = 1,
//     IsInteger = 0,
//     IsSigned = 1,
//     RequireInitialization = 1,
//     ReadCost = 1,
//     AddCost = 4,
//     MulCost = 8
//   };
// };

}  // namespace Eigen

#endif
