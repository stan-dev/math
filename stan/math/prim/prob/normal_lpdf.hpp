#ifndef STAN_MATH_PRIM_PROB_NORMAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_NORMAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {
/**
 * Return a writable column array view of the partials of a vector edge.
 * Writes through the view change `p`.
 *
 * @tparam T type of the partials (an Eigen column or row vector with
 * contiguous storage, such as an `arena_matrix`)
 * @param p partials
 * @return column array view of `p`
 */
template <typename T>
inline auto partials_column_array(T& p) {
  return Eigen::Map<Eigen::Array<value_type_t<T>, Eigen::Dynamic, 1>>(p.data(),
                                                                      p.size());
}
}  // namespace internal

/** \ingroup prob_dists
 * The log of the normal density for the specified scalar(s) given
 * the specified mean(s) and deviation(s). y, mu, or sigma can
 * each be either a scalar or a vector. Any vector inputs
 * must be the same length.
 *
 * <p>The result log probability is defined to be the sum of the
 * log probabilities for each observation/mean/deviation triple.
 *
 * @tparam T_y type of scalar
 * @tparam T_loc type of location parameter
 * @tparam T_scale type of scale parameter
 * @param y (Sequence of) scalar(s).
 * @param mu (Sequence of) location parameter(s)
 * for the normal distribution.
 * @param sigma (Sequence of) scale parameters for the normal distribution.
 * @return The log of the product of the densities.
 * @throw std::domain_error if y is NaN, the location is not finite, or the
 * scale is not positive.
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale> normal_lpdf(T_y&& y, T_loc&& mu,
                                                      T_scale&& sigma) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  static constexpr const char* function = "normal_lpdf";
  static constexpr bool y_vec_ad
      = is_autodiff_v<T_y> && !is_stan_scalar_v<std::decay_t<T_y>>;
  static constexpr bool mu_vec_ad
      = is_autodiff_v<T_loc> && !is_stan_scalar_v<std::decay_t<T_loc>>;
  static constexpr bool sigma_vec = !is_stan_scalar_v<std::decay_t<T_scale>>;
  static constexpr bool sigma_vec_ad = sigma_vec && is_autodiff_v<T_scale>;
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);
  T_y_ref y_ref = std::forward<T_y>(y);
  T_mu_ref mu_ref = std::forward<T_loc>(mu);
  T_sigma_ref sigma_ref = std::forward<T_scale>(sigma);

  // The values of autodiff containers are views, not copies.
  decltype(auto) y_val
      = to_ref_if<!is_autodiff_v<T_y>>(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref_if<!is_autodiff_v<T_loc>>(
      as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref_if<!is_autodiff_v<T_scale>>(
      as_value_column_array_or_scalar(sigma_ref));

  check_positive(function, "Scale parameter", sigma_val);

  // y and mu are checked at the early returns, and otherwise only if the sum
  // of squares is not finite: a NaN y or a non-finite mu makes it NaN or
  // infinite.
  const auto check_y_mu = [&]() {
    check_not_nan(function, "Random variable", y_val);
    check_finite(function, "Location parameter", mu_val);
  };
  if (size_zero(y_ref, mu_ref, sigma_ref)) {
    check_y_mu();
    return 0.0;
  }
  if constexpr (!include_summand<propto, T_y, T_loc, T_scale>::value) {
    check_y_mu();
    return 0.0;
  }

  auto ops_partials = make_partials_propagator(y_ref, mu_ref, sigma_ref);

  // 1 / sigma; for an autodiff vector sigma it is held in the partials of
  // sigma until they are computed.
  auto inv_sigma = [&]() {
    if constexpr (sigma_vec_ad) {
      auto inv_sigma_map
          = internal::partials_column_array(partials<2>(ops_partials));
      inv_sigma_map = inv(sigma_val);
      return inv_sigma_map;
    } else {
      return inv(sigma_val);
    }
  }();

  T_partials_return sum_sq;
  if constexpr (y_vec_ad || mu_vec_ad) {
    // The partials of mu (or of y) first hold (y - mu) / sigma.
    auto d = internal::partials_column_array(
        partials<(mu_vec_ad ? 1 : 0)>(ops_partials));
    d = (y_val - mu_val) * inv_sigma;
    sum_sq = sum(d * d);
    if (unlikely(!std::isfinite(value_of_rec(sum_sq)))) {
      check_y_mu();
    }
    // The partial of y is the negated partial of mu.
    constexpr double sign = mu_vec_ad ? 1.0 : -1.0;
    if constexpr (sigma_vec_ad) {
      for (Eigen::Index n = 0; n < d.size(); ++n) {
        const T_partials_return y_scaled = d.coeff(n);
        const T_partials_return inv_sigma_n = inv_sigma.coeff(n);
        d.coeffRef(n) = sign * (inv_sigma_n * y_scaled);
        inv_sigma.coeffRef(n)
            = inv_sigma_n * (y_scaled * y_scaled) - inv_sigma_n;
      }
    } else {
      d = sign * (inv_sigma * d);
    }
    if constexpr (mu_vec_ad && y_vec_ad) {
      internal::partials_column_array(partials<0>(ops_partials)) = -d;
    } else if constexpr (mu_vec_ad && is_autodiff_v<T_y>) {
      partials<0>(ops_partials) = -d;
    } else if constexpr (y_vec_ad && is_autodiff_v<T_loc>) {
      partials<1>(ops_partials) = -d;
    }
  } else {
    // No partial of size N other than that of sigma.
    const auto& y_scaled = (y_val - mu_val) * inv_sigma;
    sum_sq = sum(y_scaled * y_scaled);
    if (unlikely(!std::isfinite(value_of_rec(sum_sq)))) {
      check_y_mu();
    }
    if constexpr (is_autodiff_v<T_y> || is_autodiff_v<T_loc>) {
      const T_partials_return scaled_diff_sum = sum(inv_sigma * y_scaled);
      if constexpr (is_autodiff_v<T_y>) {
        partials<0>(ops_partials) = -scaled_diff_sum;
      }
      if constexpr (is_autodiff_v<T_loc>) {
        partials<1>(ops_partials) = scaled_diff_sum;
      }
    }
    if constexpr (sigma_vec_ad) {
      inv_sigma = inv_sigma * (y_scaled * y_scaled) - inv_sigma;
    }
  }

  const size_t N = max_size(y_ref, mu_ref, sigma_ref);
  if constexpr (is_autodiff_v<T_scale> && !sigma_vec) {
    partials<2>(ops_partials) = inv_sigma * (sum_sq - static_cast<double>(N));
  }

  T_partials_return logp = -0.5 * sum_sq;
  if constexpr (include_summand<propto>::value) {
    logp += NEG_LOG_SQRT_TWO_PI * N;
  }
  if constexpr (include_summand<propto, T_scale>::value) {
    logp -= sum(log(sigma_val)) * N / math::size(sigma);
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_loc, T_scale> normal_lpdf(T_y&& y, T_loc&& mu,
                                                      T_scale&& sigma) {
  return normal_lpdf<false>(std::forward<T_y>(y), std::forward<T_loc>(mu),
                            std::forward<T_scale>(sigma));
}

}  // namespace math
}  // namespace stan
#endif
