#ifndef STAN_MATH_PRIM_FUN_UB_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_UB_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/identity_constrain.hpp>
#include <stan/math/prim/fun/subtract.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the upper-bounded value for the specified unconstrained
 * matrix and upper bound.
 *
 * <p>The transform is
 *
 * <p>\f$f(x) = U - \exp(x)\f$
 *
 * <p>where \f$U\f$ is the upper bound.
 *
 * If the upper bound is positive infinity, this function
 * reduces to <code>identity_constrain(x)</code>.
 *
 * @tparam T type of Matrix
 * @tparam U type of upper bound
 * @param[in] x free Matrix.
 * @param[in] ub upper bound
 * @return matrix constrained to have upper bound
 */
template <typename T, typename L,
	  require_stan_scalar_t<L>* = nullptr>
inline promote_scalar_t<return_type_t<T, L>, T>
ub_constrain(const T& x, const L& ub) {
  if (unlikely(value_of_rec(ub) == INFTY)) {
    return x;
  }
  return subtract(ub, exp(x));
}

template <typename T, typename L,
	  require_all_matrix_t<T, L>* = nullptr,
	  require_all_not_st_var<T, L>* = nullptr>
inline auto ub_constrain(const T& x, const L& ub) {
  promote_scalar_t<return_type_t<T, L>, plain_type_t<T>> ret(x.rows(), x.cols());
  for(size_t j = 0; j < x.cols(); ++j) {
    for(size_t i = 0; i < x.rows(); ++i) {
      if (unlikely(value_of_rec(ub.coeff(i, j)) == INFTY)) {
	ret.coeffRef(i, j) = x.coeff(i, j);
      } else {
	ret.coeffRef(i, j) = ub.coeff(i, j) - exp(x.coeff(i, j));
      }
    }
  }
  return ret;
}

/**
 * Return the upper-bounded value for the specified unconstrained
 * scalar and upper bound and increment the specified log
 * probability reference with the log absolute Jacobian
 * determinant of the transform.
 *
 * <p>The transform is as specified for
 * <code>ub_constrain(T, double)</code>.  The log absolute Jacobian
 * determinant is
 *
 * <p>\f$ \log | \frac{d}{dx} -\mbox{exp}(x) + U |
 *     = \log | -\mbox{exp}(x) + 0 | = x\f$.
 *
 * If the upper bound is positive infinity, this function
 * reduces to <code>identity_constrain(x, lp)</code>.
 *
 * @tparam T type of scalar
 * @tparam U type of upper bound
 * @tparam S type of log probability
 * @param[in] x free scalar
 * @param[in] ub upper bound
 * @param[in,out] lp log density
 * @return scalar constrained to have upper bound
 */
template <typename T, typename L,
	  require_stan_scalar_t<L>* = nullptr>
inline promote_scalar_t<return_type_t<T, L>, T>
ub_constrain(const T& x, const L& ub, return_type_t<T, L>& lp) {
  if (unlikely(value_of_rec(ub) == INFTY)) {
    return x;
  }
  lp += sum(x);
  return subtract(ub, exp(x));
}

template <typename T, typename L,
	  require_all_matrix_t<T, L>* = nullptr,
	  require_all_not_st_var<T, L>* = nullptr>
inline auto ub_constrain(const T& x, const L& ub, return_type_t<T, L>& lp) {
  promote_scalar_t<return_type_t<T, L>, plain_type_t<T>> ret(x.rows(), x.cols());
  for(size_t j = 0; j < x.cols(); ++j) {
    for(size_t i = 0; i < x.rows(); ++i) {
      if (unlikely(value_of_rec(ub.coeff(i, j)) == INFTY)) {
	ret.coeffRef(i, j) = x.coeff(i, j);
      } else {
	ret.coeffRef(i, j) = ub.coeff(i, j) - exp(x.coeff(i, j));
	lp += x.coeff(i, j);
      }
    }
  }
  return ret;
}

}  // namespace math
}  // namespace stan

#endif
