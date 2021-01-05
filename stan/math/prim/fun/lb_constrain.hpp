#ifndef STAN_MATH_PRIM_FUN_LB_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_LB_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/add.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/identity_constrain.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the lower-bounded value for the specified unconstrained input
 * and specified lower bound.
 *
 * <p>The transform applied is
 *
 * <p>\f$f(x) = \exp(x) + L\f$
 *
 * <p>where \f$L\f$ is the constant lower bound.
 *
 * <p>If the lower bound is negative infinity, this function
 * reduces to <code>identity_constrain(x)</code>.
 *
 * @tparam T type of Matrix
 * @tparam L type of lower bound
 * @param[in] x Unconstrained Matrix input
 * @param[in] lb lower bound on constrained output
 * @return lower bound constrained value corresponding to inputs
 */
template <typename T, typename L, require_stan_scalar_t<L>* = nullptr>
inline promote_scalar_t<return_type_t<T, L>, T> lb_constrain(const T& x,
                                                             const L& lb) {
  if (unlikely(value_of_rec(lb) == -INFTY)) {
    return x;
  }
  return add(exp(x), lb);
}

template <typename T, typename L, require_all_matrix_t<T, L>* = nullptr,
          require_all_not_st_var<T, L>* = nullptr>
inline auto lb_constrain(const T& x, const L& lb) {
  promote_scalar_t<return_type_t<T, L>, plain_type_t<T>> ret(x.rows(),
                                                             x.cols());
  for (size_t j = 0; j < x.cols(); ++j) {
    for (size_t i = 0; i < x.rows(); ++i) {
      if (unlikely(lb.coeff(i, j) == -INFTY)) {
        ret.coeffRef(i, j) = x.coeff(i, j);
      } else {
        ret.coeffRef(i, j) = exp(x.coeff(i, j)) + lb.coeff(i, j);
      }
    }
  }
  return ret;
}

/**
 * Return the lower-bounded value for the specified unconstrained
 * input and specified lower bound, incrementing the specified
 * reference with the log absolute Jacobian determinant of the
 * transform.
 *
 * If the lower bound is negative infinity, this function
 * reduces to <code>identity_constraint(x, lp)</code>.
 *
 * @tparam T Type of Matrix
 * @tparam L type of lower bound
 * @tparam S type of log probability
 * @param[in] x unconstrained Matrix input
 * @param[in] lb lower bound on output
 * @param[in,out] lp reference to log probability to increment
 * @return lower-bound constrained value corresponding to inputs
 */
template <typename T, typename L, require_stan_scalar_t<L>* = nullptr>
inline promote_scalar_t<return_type_t<T, L>, T> lb_constrain(
    const T& x, const L& lb, return_type_t<T, L>& lp) {
  if (unlikely(value_of_rec(lb) == -INFTY)) {
    return x;
  }
  lp += sum(x);
  return add(exp(x), lb);
}

template <typename T, typename L, require_all_matrix_t<T, L>* = nullptr,
          require_all_not_st_var<T, L>* = nullptr>
inline auto lb_constrain(const T& x, const L& lb, return_type_t<T, L>& lp) {
  promote_scalar_t<return_type_t<T, L>, plain_type_t<T>> ret(x.rows(),
                                                             x.cols());
  for (size_t j = 0; j < x.cols(); ++j) {
    for (size_t i = 0; i < x.rows(); ++i) {
      if (unlikely(lb.coeff(i, j) == -INFTY)) {
        ret.coeffRef(i, j) = x.coeff(i, j);
      } else {
        ret.coeffRef(i, j) = exp(x.coeff(i, j)) + lb.coeff(i, j);
        lp += x.coeff(i, j);
      }
    }
  }
  return ret;
}

}  // namespace math
}  // namespace stan

#endif
