#ifndef STAN_MATH_PRIM_FUN_LB_FREE_HPP
#define STAN_MATH_PRIM_FUN_LB_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/identity_free.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/subtract.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the unconstrained value that produces the specified
 * lower-bound constrained value.
 *
 * @tparam T type of bounded object
 * @tparam L type of lower bound
 * @param[in] y input object
 * @param[in] lb lower bound
 * @return unconstrained value that produces the input when
 * constrained
 * @throw std::domain_error if y is lower than the lower bound
 */
template <typename T, typename L, require_not_std_vector_t<T>* = nullptr,
          require_stan_scalar_t<L>* = nullptr>
inline auto lb_free(T&& y, L&& lb) {
  if (value_of_rec(lb) == NEGATIVE_INFTY) {
    return identity_free(y, lb);
  } else {
    auto&& y_ref = to_ref(std::forward<T>(y));
    auto&& lb_ref = to_ref(std::forward<L>(lb));
    check_greater_or_equal("lb_free", "Lower bounded variable", y_ref, lb_ref);
    return eval(log(subtract(std::forward<decltype(y_ref)>(y_ref),
                             std::forward<decltype(lb_ref)>(lb_ref))));
  }
}

/**
 * Return the free matrix that corresponds to the specified
 * lower-bounded matrix with respect to the specified lower bound.
 *
 * The transform is the reverse of the `lb_constrain` transform.
 *
 * @tparam T type of matrix
 * @tparam L type of lower bound
 * @param y constrained matrix with specified lower bound
 * @param lb lower bound
 * @return unconstrained matrix with respect to lower bound
 * @throw std::invalid_argument if any element of constrained variable
 *   is greater than the lower bound.
 */
template <typename T, typename L, require_all_eigen_t<T, L>* = nullptr>
inline auto lb_free(T&& y, L&& lb) {
  auto&& y_ref = to_ref(std::forward<T>(y));
  auto&& lb_ref = to_ref(std::forward<L>(lb));
  promote_scalar_t<return_type_t<T, L>, T> ret(y.rows(), y.cols());
  for (Eigen::Index j = 0; j < y.cols(); ++j) {
    for (Eigen::Index i = 0; i < y.rows(); ++i) {
      ret.coeffRef(i, j) = lb_free(y_ref.coeff(i, j), lb_ref.coeff(i, j));
    }
  }
  return ret;
}

/**
 * Return the free variable that corresponds to the specified
 * lower-bounded variable with respect to the specified lower bound.
 *
 * The transform is the reverse of the `lb_constrain` transform.
 *
 * @tparam T type of constrained variable
 * @tparam L type of lower bound
 * @param y constrained variable with specified lower bound
 * @param lb lower bound
 * @return unconstrained variable with respect to lower bound
 * @throw std::invalid_argument if any element of constrained variable
 *   is greater than the lower bound.
 */
template <typename T, typename L, require_not_std_vector_t<L>* = nullptr>
inline auto lb_free(const std::vector<T> y, const L& lb) {
  auto&& lb_ref = to_ref(lb);
  std::vector<decltype(lb_free(y[0], lb))> ret(y.size());
  std::transform(y.begin(), y.end(), ret.begin(),
                 [&lb](auto&& yy) { return lb_free(yy, lb); });
  return ret;
}

/**
 * Return the free variable that corresponds to the specified
 * lower-bounded variable with respect to the specified lower bound.
 *
 * The transform is the reverse of the `lb_constrain` transform.
 *
 * @tparam T type of constrained variable
 * @tparam L type of lower bound
 * @param y constrained variable with specified lower bound
 * @param lb lower bound
 * @return unconstrained variable with respect to lower bound
 * @throw std::invalid_argument if any element of constrained variable
 *   is greater than the lower bound.
 */
template <typename T, typename L>
inline auto lb_free(const std::vector<T> y, const std::vector<L>& lb) {
  std::vector<decltype(lb_free(y[0], lb[0]))> ret(y.size());
  std::transform(y.begin(), y.end(), lb.begin(), ret.begin(),
                 [](auto&& yy, auto&& lbb) { return lb_free(yy, lbb); });
  return ret;
}

}  // namespace math
}  // namespace stan
#endif
