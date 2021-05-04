#ifndef STAN_MATH_PRIM_FUN_UB_FREE_HPP
#define STAN_MATH_PRIM_FUN_UB_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/identity_free.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the free scalar that corresponds to the specified
 * upper-bounded value with respect to the specified upper bound.
 *
 * <p>The transform is the reverse of the
 * <code>ub_constrain(T, double)</code> transform,
 *
 * <p>\f$f^{-1}(y) = \log -(y - U)\f$
 *
 * <p>where \f$U\f$ is the upper bound.
 *
 * @tparam T type of scalar
 * @tparam U type of upper bound
 * @param y constrained scalar with specified upper bound
 * @param ub upper bound
 * @return unconstrained scalar with respect to upper bound
 * @throw std::invalid_argument if constrained scalar is greater
 *   than the upper bound.
 */
template <typename T, typename U, require_not_std_vector_t<T>* = nullptr,
          require_stan_scalar_t<U>* = nullptr>
inline auto ub_free(T&& y, U&& ub) {
  if (value_of_rec(ub) == INFTY) {
    return identity_free(y, ub);
  } else {
    auto&& y_ref = to_ref(std::forward<T>(y));
    auto&& ub_ref = to_ref(std::forward<U>(ub));
    check_less_or_equal("ub_free", "Upper bounded variable", value_of(y_ref),
                        value_of(ub_ref));
    return eval(log(subtract(std::forward<decltype(ub_ref)>(ub_ref),
                             std::forward<decltype(y_ref)>(y_ref))));
  }
}

/**
 * Return the free matrix that corresponds to the specified
 * upper-bounded matrix with respect to the specified upper bound.
 *
 * The transform is the reverse of the `ub_constrain` transform.
 *
 * @tparam T type of matrix
 * @tparam U type of upper bound
 * @param y constrained matrix with specified upper bound
 * @param ub upper bound
 * @return unconstrained matrix with respect to upper bound
 * @throw std::invalid_argument if any element of constrained variable
 *   is greater than the upper bound.
 */
template <typename T, typename U, require_all_eigen_t<T, U>* = nullptr>
inline auto ub_free(T&& y, U&& ub) {
  check_matching_dims("ub_free", "y", y, "ub", ub);
  auto&& y_ref = to_ref(std::forward<T>(y));
  auto&& ub_ref = to_ref(std::forward<U>(ub));
  promote_scalar_t<return_type_t<T, U>, T> ret(y.rows(), y.cols());
  for (Eigen::Index j = 0; j < y.cols(); ++j) {
    for (Eigen::Index i = 0; i < y.rows(); ++i) {
      ret.coeffRef(i, j) = ub_free(y_ref.coeff(i, j), ub_ref.coeff(i, j));
    }
  }
  return ret;
}

/**
 * Return the free variable that corresponds to the specified
 * upper-bounded variable with respect to the specified upper bound.
 *
 * The transform is the reverse of the `ub_constrain` transform.
 *
 * @tparam T type of constrained variable
 * @tparam U type of upper bound
 * @param y constrained variable with specified upper bound
 * @param ub upper bound
 * @return unconstrained variable with respect to upper bound
 * @throw std::invalid_argument if any element of constrained variable
 *   is greater than the upper bound.
 */
template <typename T, typename U, require_not_std_vector_t<U>* = nullptr>
inline auto ub_free(const std::vector<T> y, const U& ub) {
  auto&& ub_ref = to_ref(ub);
  std::vector<decltype(ub_free(y[0], ub))> ret(y.size());
  for (Eigen::Index i = 0; i < y.size(); ++i) {
    ret[i] = ub_free(y[i], ub_ref);
  }
  return ret;
}

/**
 * Return the free variable that corresponds to the specified
 * upper-bounded variable with respect to the specified upper bound.
 *
 * The transform is the reverse of the `ub_constrain` transform.
 *
 * @tparam T type of constrained variable
 * @tparam U type of upper bound
 * @param y constrained variable with specified upper bound
 * @param ub upper bound
 * @return unconstrained variable with respect to upper bound
 * @throw std::invalid_argument if any element of constrained variable
 *   is greater than the upper bound.
 */
template <typename T, typename U>
inline auto ub_free(const std::vector<T> y, const std::vector<U>& ub) {
  check_matching_dims("ub_free", "y", y, "ub", ub);
  std::vector<decltype(ub_free(y[0], ub[0]))> ret(y.size());
  for (Eigen::Index i = 0; i < y.size(); ++i) {
    ret[i] = ub_free(y[i], ub[i]);
  }
  return ret;
}
}  // namespace math
}  // namespace stan
#endif
