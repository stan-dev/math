#ifndef STAN_MATH_PRIM_FUN_LB_FREE_HPP
#define STAN_MATH_PRIM_FUN_LB_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/identity_free.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the unconstrained value that produces the specified
 * lower-bound constrained value.
 *
 * If the lower bound is negative infinity, it is ignored and
 * the function reduces to <code>identity_free(x)</code>.
 *
 * @tparam T type of scalar
 * @tparam U type of lower bound
 * @param[in] x input scalar
 * @param[in] lb lower bound
 * @return unconstrained value that produces the input when
 * constrained
 * @throw std::domain_error if x is lower than the lower bound
 */
template <typename T, typename U, require_stan_scalar_t<T>* = nullptr>
inline auto lb_free(T&& x, U&& lb) {
  using std::log;
  if (lb == NEGATIVE_INFTY) {
    return identity_free(x, lb);
  }
  check_greater_or_equal("lb_free", "Lower bounded variable", x, lb);
  return log(x - lb);
}

/**
 * Return the unconstrained value that produces the specified
 * lower-bound constrained value.
 *
 * If the lower bound is negative infinity, it is ignored and
 * the function reduces to <code>identity_free(x)</code>.
 *
 * @tparam EigT type derived from `EigenBase`
 * @tparam U type of lower bound
 * @param[in] x input scalar
 * @param[in] lb lower bound
 * @return unconstrained value that produces the input when
 * constrained
 * @throw std::domain_error if x is lower than the lower bound
 */
template <typename EigT, typename U, require_eigen_t<EigT>* = nullptr>
inline auto lb_free(EigT&& x, U&& lb) {
  return x.unaryExpr([&lb](auto&& y_iter){
    return lb_free(y_iter, lb);
  }).eval();
}

/**
 * Return the unconstrained value that produces the specified
 * lower-bound constrained value.
 *
 * If the lower bound is negative infinity, it is ignored and
 * the function reduces to <code>identity_free(x)</code>.
 *
 * @tparam Vec type of standard vector
 * @tparam U type of lower bound
 * @param[in] x input scalar
 * @param[in] lb lower bound
 * @return unconstrained value that produces the input when
 * constrained
 * @throw std::domain_error if x is lower than the lower bound
 */
template <typename Vec, typename U, require_std_vector_t<Vec>* = nullptr>
inline auto lb_free(Vec&& x, U&& lb) {
  std::vector<return_type_t<Vec, U>> ret_x(x.size());
  std::transform(x.begin(), x.end(), ret_x.begin(), [&lb](auto&& x_iter) {
    return lb_free(x_iter, lb);
  });
  return ret_x;
}

}  // namespace math
}  // namespace stan
#endif
