#ifndef STAN_MATH_FWD_FUN_EIGEN_NUMTRAITS_HPP
#define STAN_MATH_FWD_FUN_EIGEN_NUMTRAITS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/core.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/core/std_numeric_limits.hpp>
#include <limits>

namespace Eigen {

/**
 * Numerical traits template override for Eigen for automatic
 * gradient variables.
 */
template <typename T>
struct NumTraits<stan::math::fvar<T>> : GenericNumTraits<stan::math::fvar<T>> {
  enum {
    /**
     * stan::math::fvar requires initialization
     */
    RequireInitialization = 1,

    /**
     * twice the cost to copy a double
     */
    ReadCost = 2 * NumTraits<double>::ReadCost,

    /**
     * 2 * AddCost
     * <br>
     * (a + b) = a + b
     * <br>
     * (a + b)' = a' + b'
     */
    AddCost = 2 * NumTraits<T>::AddCost,

    /**
     * 3 * MulCost + AddCost
     * <br>
     * (a * b) = a * b
     * <br>
     * (a * b)' = a' * b + a * b'
     */
    MulCost = 3 * NumTraits<T>::MulCost + NumTraits<T>::AddCost
  };

  /**
   * Return the number of decimal digits that can be represented
   * without change.  Delegates to
   * <code>std::numeric_limits<double>::digits10()</code>.
   */
  static int digits10() { return std::numeric_limits<double>::digits10; }
};


#define STAN_EIGEN_SCALAR_BINARY_OP_TRAITS(lhs, rhs, returntype) \
 /*! Traits specialization for Eigen binary operations for autodiff lhs */ \
 /*! rhs arguments. */ \
 /*! @tparam T value and tangent type of autodiff variable */ \
 /*! @tparam BinaryOp type of binary operation for which traits are */ \
 /*! defined */ \
template <typename T, typename BinaryOp> \
struct ScalarBinaryOpTraits<lhs, rhs, BinaryOp> { \
  using ReturnType = returntype; \
}; \
/*! Traits specialization for Eigen binary operations for autodiff rhs */ \
/*! lhs arguments. */ \
/*! @tparam T value and tangent type of autodiff variable */ \
/*! @tparam BinaryOp type of binary operation for which traits are */ \
/*! defined */ \
template <typename T, typename BinaryOp> \
struct ScalarBinaryOpTraits<rhs, lhs, BinaryOp> { \
  using ReturnType = returntype; \
};

STAN_EIGEN_SCALAR_BINARY_OP_TRAITS(stan::math::fvar<T>, double, stan::math::fvar<T>);
STAN_EIGEN_SCALAR_BINARY_OP_TRAITS(stan::math::fvar<T>, std::complex<double>, std::complex<stan::math::fvar<T>>);
STAN_EIGEN_SCALAR_BINARY_OP_TRAITS(std::complex<stan::math::fvar<T>>, double, std::complex<stan::math::fvar<T>>);
STAN_EIGEN_SCALAR_BINARY_OP_TRAITS(std::complex<double>, std::complex<stan::math::fvar<T>>, std::complex<stan::math::fvar<T>>);

#undef STAN_EIGEN_SCALAR_BINARY_OP_TRAITS

template <typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<stan::math::fvar<T>, stan::math::fvar<T>, BinaryOp> {
  using ReturnType = stan::math::fvar<T>;
};

template <typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<std::complex<stan::math::fvar<T>>, std::complex<stan::math::fvar<T>>, BinaryOp> {
  using ReturnType = std::complex<stan::math::fvar<T>>;
};

}  // namespace Eigen
#endif
