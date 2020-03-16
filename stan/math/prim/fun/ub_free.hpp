#ifndef STAN_MATH_PRIM_FUN_UB_FREE_HPP
#define STAN_MATH_PRIM_FUN_UB_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/identity_free.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <algorithm>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the free scalar that corresponds to the specified
 * upper-bounded value with respect to the specified upper bound.
 *
 * <p>The transform is the reverse of the
 * <code>ub_constrain(T, double)</code> transform,
 *
 * <p>\f$f^{-1}(x) = \log -(x - U)\f$
 *
 * <p>where \f$U\f$ is the upper bound.
 *
 * If the upper bound is positive infinity, this function
 * reduces to <code>identity_free(x, ub)</code>.
 *
 * @tparam T type of scalar
 * @tparam U type of upper bound
 * @param x constrained scalar with specified upper bound
 * @param ub upper bound
 * @return unconstrained scalar with respect to upper bound
 * @throw std::invalid_argument if constrained scalar is greater
 *   than the upper bound.
 */
template <typename T, typename U, require_stan_scalar_t<T>* = nullptr>
inline auto ub_free(T&& x, U&& ub) {
  using std::log;
  if (ub == INFTY) {
    return identity_free(x, ub);
  }
  check_less_or_equal("ub_free", "Upper bounded variable", x, ub);
  return log(ub - x);
}

/**
 * Return the free scalar that corresponds to the specified
 * upper-bounded value with respect to the specified upper bound.
 *
 * <p>The transform is the reverse of the
 * <code>ub_constrain(T, double)</code> transform,
 *
 * <p>\f$f^{-1}(x) = \log -(x - U)\f$
 *
 * <p>where \f$U\f$ is the upper bound.
 *
 * If the upper bound is positive infinity, this function
 * reduces to <code>identity_free(x, ub)</code>.
 *
 * @tparam EigT Type derived from EigenBase
 * @tparam U type of upper bound
 * @param x constrained scalar with specified upper bound
 * @param ub upper bound
 * @return unconstrained scalar with respect to upper bound
 * @throw std::invalid_argument if constrained scalar is greater
 *   than the upper bound.
 */
template <typename EigT, typename U, require_eigen_t<EigT>* = nullptr>
inline auto ub_free(EigT&& x, U&& ub) {
  using std::log;
  if (ub == INFTY) {
    return identity_free(std::forward<EigT>(x), ub);
  }
  check_less_or_equal("ub_free", "Upper bounded variable", x, ub);
  return (ub - x.array()).log().matrix().eval();
}

/**
 * Return the free scalar that corresponds to the specified
 * upper-bounded value with respect to the specified upper bound.
 *
 * <p>The transform is the reverse of the
 * <code>ub_constrain(T, double)</code> transform,
 *
 * <p>\f$f^{-1}(x) = \log -(x - U)\f$
 *
 * <p>where \f$U\f$ is the upper bound.
 *
 * If the upper bound is positive infinity, this function
 * reduces to <code>identity_free(x)</code>.
 *
 * @tparam Vec type of standard vector
 * @tparam U type of upper bound
 * @param x constrained scalar with specified upper bound
 * @param ub upper bound
 * @return unconstrained scalar with respect to upper bound
 * @throw std::invalid_argument if constrained scalar is greater
 *   than the upper bound.
 */
template <typename Vec, typename U, require_std_vector_t<Vec>* = nullptr>
inline auto ub_free(Vec&& x, U&& ub) {
  if (ub == INFTY) {
    return identity_free(std::forward<Vec>(x), ub);
  }
  check_less_or_equal("ub_free", "Upper bounded variable", x, ub);
  std::vector<return_type_t<Vec, U>> ret_x(x.size());
  std::transform(x.begin(), x.end(), ret_x.begin(),
                 [&ub](auto&& x_iter) { return log(ub - x_iter); });
  return ret_x;
}

}  // namespace math
}  // namespace stan
#endif
