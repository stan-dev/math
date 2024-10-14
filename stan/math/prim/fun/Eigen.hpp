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
#include <stan/math/prim/core/complex_base.hpp>
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
  struct ScalarBinaryOpTraits<int, stan::math::complex<double>, BinaryOp> {
    using ReturnType = stan::math::complex<double>;
  };

  /**
   * Traits specialization for Eigen binary operations for complex
   * `double` and `int` arguments.
   *
   * @tparam BinaryOp type of binary operation for which traits are
   * defined
   */
  template <typename BinaryOp>
  struct ScalarBinaryOpTraits<stan::math::complex<double>, int, BinaryOp> {
    using ReturnType = stan::math::complex<double>;
  };

  template <typename Real_>
struct NumTraits<stan::math::complex<Real_> > : GenericNumTraits<stan::math::complex<Real_> > {
  typedef Real_ Real;
  typedef typename NumTraits<Real_>::Literal Literal;
  enum {
    IsComplex = 1,
    IsSigned = NumTraits<Real_>::IsSigned,
    RequireInitialization = NumTraits<Real_>::RequireInitialization,
    ReadCost = 2 * NumTraits<Real_>::ReadCost,
    AddCost = 2 * NumTraits<Real>::AddCost,
    MulCost = 4 * NumTraits<Real>::MulCost + 2 * NumTraits<Real>::AddCost
  };

  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR static inline Real epsilon() { return NumTraits<Real>::epsilon(); }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR static inline Real dummy_precision() { return NumTraits<Real>::dummy_precision(); }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR static inline int digits10() { return NumTraits<Real>::digits10(); }
  EIGEN_DEVICE_FUNC EIGEN_CONSTEXPR static inline int max_digits10() { return NumTraits<Real>::max_digits10(); }
};



}  // namespace Eigen

#endif
