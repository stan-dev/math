#ifndef STAN_MATH_PRIM_FUNCTOR_APPLY_SCALAR_UNARY_HPP
#define STAN_MATH_PRIM_FUNCTOR_APPLY_SCALAR_UNARY_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_complex.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/is_vector_like.hpp>
#include <stan/math/prim/meta/plain_type.hpp>
#include <utility>
#include <vector>

namespace stan {
namespace math {

/**
 * Base template class for vectorization of unary scalar functions
 * defined by a template class <code>F</code> to a scalar,
 * standard library vector, or Eigen dense matrix expression
 * template.
 *
 * <p>Specializations of this base define applications to scalars
 * (primitive or autodiff in the corresponding autodiff library
 * directories) or to standard library vectors of vectorizable
 * types (primitives, Eigen dense matrix expressions, or further
 * standard vectors).
 *
 * <p>Each specialization must define the typedef
 * <code>return_t</code> for the vectorized return type and the
 * function <code>apply</code> which defines the vectorization or
 * base application of the function defined statically by the
 * class F.  The function definition class F defines a static
 * function <code>fun()</code>, which defines the function's
 * behavior on scalars.
 *
 * @tparam F Type of function to apply.
 * @tparam T Type of argument to which function is applied.
 */
template <typename F, typename T, typename Enable = void>
struct apply_scalar_unary;

/**
 *
 * Template specialization for vectorized functions applying to
 * Eigen matrix arguments.
 *
 * @tparam F Type of function to apply.
 * @tparam T Type of argument to which function is applied.
 */
template <typename F, typename T>
struct apply_scalar_unary<F, T, require_eigen_t<T>> {
  /**
   * Return the result of applying the function defined by the
   * template parameter F to the specified matrix argument.
   *
   * @param x Matrix to which operation is applied.
   * @return Componentwise application of the function specified
   * by F to the specified matrix.
   */
  static inline auto apply(const T& x) {
    return x.unaryExpr([](auto&& x) {
      return apply_scalar_unary<F, std::decay_t<decltype(x)>>::apply(x);
    });
  }

  /**
   * Return type for applying the function elementwise to a matrix
   * expression template of type T.
   */
  using return_t = std::decay_t<decltype(
      apply_scalar_unary<F, T>::apply(std::declval<T>()))>;
};

/**
 * Template specialization for vectorized functions applying to
 * double arguments.
 *
 * @tparam F Type of function defining static apply function.
 */
template <typename F, typename T>
struct apply_scalar_unary<F, T, require_floating_point_t<T>> {
  /**
   * The return type, double.
   */
  using return_t = std::decay_t<decltype(F::fun(std::declval<T>()))>;

  /**
   * Apply the function specified by F to the specified argument.
   * This is defined through a direct application of
   * <code>F::fun()</code>, which must be defined for double
   * arguments.
   *
   * @param x Argument scalar.
   * @return Result of applying F to the scalar.
   */
  static inline auto apply(T x) { return F::fun(x); }
};

/**
 * Template specialization for vectorized functions applying to
 * complex arguments.
 *
 * @tparam F Type of function defining static apply function.
 */
template <typename F, typename T>
struct apply_scalar_unary<F, T, require_complex_t<T>> {
  /**
   * Apply the function specified by F to the specified argument.
   * This is defined through a direct application of
   * <code>F::fun()</code>, which must be defined for double
   * arguments.
   *
   * @param x Argument scalar.
   * @return Result of applying F to the scalar.
   */
  static inline auto apply(const T& x) { return F::fun(x); }
  /**
   * The return type
   */
  using return_t = std::decay_t<decltype(F::fun(std::declval<T>()))>;
};

/**
 * Template specialization for vectorized functions applying to
 * integer arguments.  Although the argument is integer, the
 * return type is specified as double.  This allows promotion of
 * integers to doubles in vectorized functions, or in containers.
 *
 * @tparam F Type of function defining static apply function.
 */
template <typename F, typename T>
struct apply_scalar_unary<F, T, require_integral_t<T>> {
  /**
   * Apply the function specified by F to the specified argument.
   * This is defined through a direct application of
   * <code>F::fun()</code>, which must be defined for double
   * arguments.
   *
   * @param x Argument scalar.
   * @return Result of applying F to the scalar.
   */
  static inline auto apply(T x) { return F::fun(x); }
  /**
   * The return type, double.
   */
  using return_t = std::decay_t<decltype(F::fun(std::declval<double>()))>;
};

/**
 * Template specialization for vectorized functions applying to
 * standard vector containers.  The lowest-level scalar type of
 * the argument will determine the return type.  Integers are
 * promoted to double values.
 *
 * @tparam F Class defining a static apply function.
 * @tparam T Type of element contained in standard vector.
 */
template <typename F, typename T>
struct apply_scalar_unary<F, std::vector<T>> {
  /**
   * Return type, which is calculated recursively as a standard
   * vector of the return type of the contained type T.
   */
  using return_t = typename std::vector<
      plain_type_t<typename apply_scalar_unary<F, T>::return_t>>;

  /**
   * Apply the function specified by F elementwise to the
   * specified argument.  This is defined recursively through this
   * class applied to elements of type T.
   *
   * @param x Argument container.
   * @return Elementwise application of F to the elements of the
   * container.
   */
  static inline auto apply(const std::vector<T>& x) {
    return_t fx(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
      fx[i] = apply_scalar_unary<F, T>::apply(x[i]);
    }
    return fx;
  }
};

}  // namespace math
}  // namespace stan
#endif
