#ifndef STAN_MATH_FWD_FUN_EIGEN_NUMTRAITS_HPP
#define STAN_MATH_FWD_FUN_EIGEN_NUMTRAITS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/core.hpp>
#include <stan/math/fwd/fun/read_fvar.hpp>
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

/**
 * Traits specialization for Eigen binary operations for autodiff and
 * `double` arguments.
 *
 * @tparam T value and tangent type of autodiff variable
 * @tparam BinaryOp type of binary operation for which traits are
 * defined
 */
template <typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<stan::math::fvar<T>, double, BinaryOp> {
  using ReturnType = stan::math::fvar<T>;
};

/**
 * Traits specialization for Eigen binary operations for `double` and
 * autodiff arguments.
 *
 * @tparam T value and tangent type of autodiff variable
 * @tparam BinaryOp type of binary operation for which traits are
 * defined
 */
template <typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<double, stan::math::fvar<T>, BinaryOp> {
  using ReturnType = stan::math::fvar<T>;
};

/**
 * Traits specialization for Eigen binary operations for autodiff and
 * `int` arguments.
 *
 * @tparam T value and tangent type of autodiff variable
 * @tparam BinaryOp type of binary operation for which traits are
 * defined
 */
template <typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<stan::math::fvar<T>, int, BinaryOp> {
  using ReturnType = stan::math::fvar<T>;
};

/**
 * Traits specialization for Eigen binary operations for `int` and
 * autodiff arguments.
 *
 * @tparam T value and tangent type of autodiff variable
 * @tparam BinaryOp type of binary operation for which traits are
 * defined
 */
template <typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<int, stan::math::fvar<T>, BinaryOp> {
  using ReturnType = stan::math::fvar<T>;
};

/**
 * Traits specialization for Eigen binary operations for `double` and
 * complex autodiff arguments.
 *
 * @tparam T value and tangent type of autodiff variable
 * @tparam BinaryOp type of binary operation for which traits are
 * defined
 */
template <typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<double, std::complex<stan::math::fvar<T>>,
                            BinaryOp> {
  using ReturnType = std::complex<stan::math::fvar<T>>;
};

/**
 * Traits specialization for Eigen binary operations for complex
 * autodiff and `double` arguments.
 *
 * @tparam T value and tangent type of autodiff variable
 * @tparam BinaryOp type of binary operation for which traits are
 * defined
 */
template <typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<std::complex<stan::math::fvar<T>>, double,
                            BinaryOp> {
  using ReturnType = std::complex<stan::math::fvar<T>>;
};

/**
 * Traits specialization for Eigen binary operations for `int` and
 * complex autodiff arguments.
 *
 * @tparam T value and tangent type of autodiff variable
 * @tparam BinaryOp type of binary operation for which traits are
 * defined
 */
template <typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<int, std::complex<stan::math::fvar<T>>, BinaryOp> {
  using ReturnType = std::complex<stan::math::fvar<T>>;
};

/**
 * Traits specialization for Eigen binary operations for complex
 * autodiff and `int` arguments.
 *
 * @tparam T value and tangent type of autodiff variable
 * @tparam BinaryOp type of binary operation for which traits are
 * defined
 */
template <typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<std::complex<stan::math::fvar<T>>, int, BinaryOp> {
  using ReturnType = std::complex<stan::math::fvar<T>>;
};

/**
 * Traits specialization for Eigen binary operations for autodiff
 * and complex `double` arguments.
 *
 * @tparam T value and tangent type of autodiff variable
 * @tparam BinaryOp type of binary operation for which traits are
 * defined
 */
template <typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<stan::math::fvar<T>, std::complex<double>,
                            BinaryOp> {
  using ReturnType = std::complex<stan::math::fvar<T>>;
};

/**
 * Traits specialization for Eigen binary operations for complex
 * `double` and autodiff arguments.
 *
 * @tparam T value and tangent type of autodiff variable
 * @tparam BinaryOp type of binary operation for which traits are
 * defined
 */
template <typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<std::complex<double>, stan::math::fvar<T>,
                            BinaryOp> {
  using ReturnType = std::complex<stan::math::fvar<T>>;
};

/**
 * Traits specialization for Eigen binary operations for complex
 * `double` and complex autodiff arguments.
 *
 * @tparam T value and tangent type of autodiff variable
 * @tparam BinaryOp type of binary operation for which traits are
 * defined
 */
template <typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<std::complex<double>,
                            std::complex<stan::math::fvar<T>>, BinaryOp> {
  using ReturnType = std::complex<stan::math::fvar<T>>;
};

/**
 * Traits specialization for Eigen binary operations for complex
 * autodiff and complex `double` arguments.
 *
 * @tparam T value and tangent type of autodiff variable
 * @tparam BinaryOp type of binary operation for which traits are
 * defined
 */
template <typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<std::complex<stan::math::fvar<T>>,
                            std::complex<double>, BinaryOp> {
  using ReturnType = std::complex<stan::math::fvar<T>>;
};

namespace internal {

/**
 * Enable linear access of inputs when using read_fvar.
 */
template <typename EigFvar, typename EigOut>
struct functor_has_linear_access<
    stan::math::read_fvar_functor<EigFvar, EigOut>> {
  enum { ret = 1 };
};

}  // namespace internal
}  // namespace Eigen
#endif
