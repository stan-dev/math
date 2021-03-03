#ifndef STAN_MATH_PRIM_FUN_OFFSET_MULTIPLIER_FREE_HPP
#define STAN_MATH_PRIM_FUN_OFFSET_MULTIPLIER_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/identity_free.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the unconstrained variable that transforms to the
 * specified offset and multiplier constrained variable given the specified
 * offset and multiplier.
 *
 * <p>The transform in <code>locmultiplier_constrain(T, double, double)</code>,
 * is reversed by the reverse affine transformation,
 *
 * <p>\f$f^{-1}(y) = \frac{y - L}{S}\f$
 *
 * where \f$L\f$ and \f$S\f$ are the offset and multiplier.
 *
 * <p>If the offset is zero and multiplier is one,
 * this function reduces to  <code>identity_free(y)</code>.
 *
 * @tparam T type of constrained variable
 * @tparam L type of offset
 * @tparam S type of multiplier
 * @param y constrained value
 * @param[in] mu offset of constrained output
 * @param[in] sigma multiplier of constrained output
 * @return the unconstrained variable that transforms to the given constrained
 *  variable given the offset and multiplier
 * @throw std::domain_error if sigma <= 0
 * @throw std::domain_error if mu is not finite
 * @throw std::invalid_argument if non-scalar arguments don't match in size
 */
template <typename T, typename M, typename S>
inline auto offset_multiplier_free(const T& y, const M& mu, const S& sigma) {
  auto&& mu_ref = to_ref(mu);
  auto&& sigma_ref = to_ref(sigma);
  if (is_matrix<T>::value && is_matrix<M>::value) {
    check_matching_dims("offset_multiplier_constrain", "y", y, "mu", mu);
  }
  if (is_matrix<T>::value && is_matrix<S>::value) {
    check_matching_dims("offset_multiplier_constrain", "y", y, "sigma", sigma);
  } else if (is_matrix<M>::value && is_matrix<S>::value) {
    check_matching_dims("offset_multiplier_constrain", "mu", mu, "sigma",
                        sigma);
  }

  check_finite("offset_multiplier_constrain", "offset", value_of(mu_ref));
  check_positive_finite("offset_multiplier_constrain", "multiplier",
                        value_of(sigma_ref));
  return divide(subtract(y, mu_ref), sigma_ref);
}

/**
 * Overload for array of x and non-array mu and sigma
 */
template <typename T, typename M, typename S,
          require_all_not_std_vector_t<M, S>* = nullptr>
inline auto offset_multiplier_free(const std::vector<T>& x, const M& mu,
                                   const S& sigma) {
  std::vector<plain_type_t<decltype(offset_multiplier_free(x[0], mu, sigma))>>
      ret;
  ret.reserve(x.size());
  const auto& mu_ref = to_ref(mu);
  const auto& sigma_ref = to_ref(sigma);
  for (size_t i = 0; i < x.size(); ++i) {
    ret.emplace_back(offset_multiplier_free(x[i], mu_ref, sigma_ref));
  }
  return ret;
}

/**
 * Overload for array of x and sigma and non-array mu
 */
template <typename T, typename M, typename S,
          require_not_std_vector_t<M>* = nullptr>
inline auto offset_multiplier_free(const std::vector<T>& x, const M& mu,
                                   const std::vector<S>& sigma) {
  check_matching_dims("offset_multiplier_free", "x", x, "sigma", sigma);
  std::vector<
      plain_type_t<decltype(offset_multiplier_free(x[0], mu, sigma[0]))>>
      ret;
  ret.reserve(x.size());
  const auto& mu_ref = to_ref(mu);
  for (size_t i = 0; i < x.size(); ++i) {
    ret.emplace_back(offset_multiplier_free(x[i], mu_ref, sigma[i]));
  }
  return ret;
}

/**
 * Overload for array of x and mu and non-array sigma
 */
template <typename T, typename M, typename S,
          require_not_std_vector_t<S>* = nullptr>
inline auto offset_multiplier_free(const std::vector<T>& x,
                                   const std::vector<M>& mu, const S& sigma) {
  check_matching_dims("offset_multiplier_free", "x", x, "mu", mu);
  std::vector<
      plain_type_t<decltype(offset_multiplier_free(x[0], mu[0], sigma))>>
      ret;
  ret.reserve(x.size());
  const auto& sigma_ref = to_ref(sigma);
  for (size_t i = 0; i < x.size(); ++i) {
    ret.emplace_back(offset_multiplier_free(x[i], mu[i], sigma_ref));
  }
  return ret;
}

/**
 * Overload for array of x, mu, and sigma
 */
template <typename T, typename M, typename S>
inline auto offset_multiplier_free(const std::vector<T>& x,
                                   const std::vector<M>& mu,
                                   const std::vector<S>& sigma) {
  check_matching_dims("offset_multiplier_free", "x", x, "mu", mu);
  check_matching_dims("offset_multiplier_free", "x", x, "sigma", sigma);
  std::vector<
      plain_type_t<decltype(offset_multiplier_free(x[0], mu[0], sigma[0]))>>
      ret;
  ret.reserve(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret.emplace_back(offset_multiplier_free(x[i], mu[i], sigma[i]));
  }
  return ret;
}

}  // namespace math
}  // namespace stan
#endif
