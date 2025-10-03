#ifndef STAN_MATH_PRIM_CONSTRAINT_POSITIVE_CONSTRAIN_HPP
#define STAN_MATH_PRIM_CONSTRAINT_POSITIVE_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the positive value for the specified unconstrained input.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = \exp(x)\f$.
 *
 * @param x Arbitrary input scalar or container.
 * @return Input transformed to be positive.
 */
template <typename T>
inline auto positive_constrain(T&& x) {
  return exp(std::forward<T>(x));
}

/**
 * Return the positive value for the specified unconstrained input,
 * incrementing the scalar reference with the log absolute
 * Jacobian determinant.
 *
 * <p>See <code>positive_constrain(T)</code> for details
 * of the transform.  The log absolute Jacobian determinant is
 *
 * <p>\f$\log | \frac{d}{dx} \mbox{exp}(x) |
 *    = \log | \mbox{exp}(x) | =  x\f$.
 *
 * @tparam T type of unconstrained value
 * @param x unconstrained value or container
 * @param lp log density reference.
 * @return positive constrained version of unconstrained value(s)
 */
template <typename T, typename S>
inline auto positive_constrain(T&& x, S& lp) {
  auto&& x_ref = to_ref(std::forward<T>(x));
  lp += sum(x_ref);
  return exp(std::forward<decltype(x_ref)>(x_ref));
}

/**
 * Return the positive value for the specified unconstrained input. If the
 * `Jacobian` parameter is `true`, the log density accumulator is incremented
 * with the log absolute Jacobian determinant of the transform.  All of the
 * transforms are specified with their Jacobians in the *Stan Reference Manual*
 * chapter Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A type inheriting from `Eigen::EigenBase`, a `var_value` with inner
 * type inheriting from `Eigen::EigenBase`, a standard vector, or a scalar
 * @tparam Lp A scalar type for the lp argument. The scalar type of T should be
 * convertable to this.
 * @param x unconstrained value or container
 * @param[in, out] lp log density accumulator
 * @return positive constrained version of unconstrained value(s)
 */
template <bool Jacobian, typename T, typename Lp,
          require_not_std_vector_t<T>* = nullptr,
          require_convertible_t<return_type_t<T>, Lp>* = nullptr>
inline auto positive_constrain(T&& x, Lp& lp) {
  if constexpr (Jacobian) {
    return positive_constrain(std::forward<T>(x), lp);
  } else {
    return positive_constrain(std::forward<T>(x));
  }
}

/**
 * Return the positive value for the specified unconstrained input. If the
 * `Jacobian` parameter is `true`, the log density accumulator is incremented
 * with the log absolute Jacobian determinant of the transform.  All of the
 * transforms are specified with their Jacobians in the *Stan Reference Manual*
 * chapter Constraint Transforms.
 *
 * @tparam Jacobian if `true`, increment log density accumulator with log
 * absolute Jacobian determinant of constraining transform
 * @tparam T A standard vector with inner type inheriting from
 * `Eigen::EigenBase`, a `var_value` with inner type inheriting from
 * `Eigen::EigenBase`, a standard vector, or a scalar
 * @tparam Lp A scalar type for the lp argument. The scalar type of T should be
 * convertable to this.
 * @param x unconstrained value or container
 * @param[in, out] lp log density accumulator
 * @return positive constrained version of unconstrained value(s)
 */
template <bool Jacobian, typename T, typename Lp,
          require_std_vector_t<T>* = nullptr,
          require_convertible_t<return_type_t<T>, Lp>* = nullptr>
inline auto positive_constrain(T&& x, Lp& lp) {
  return apply_vector_unary<T>::apply(std::forward<T>(x), [&lp](auto&& v) {
    return positive_constrain<Jacobian>(std::forward<decltype(v)>(v), lp);
  });
}

}  // namespace math
}  // namespace stan

#endif
