#ifndef STAN_MATH_PRIM_MAT_PROB_STD_BINORMAL_LCDF_HPP
#define STAN_MATH_PRIM_MAT_PROB_STD_BINORMAL_LCDF_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/mat/meta/operands_and_partials.hpp>
#include <stan/math/prim/mat/meta/vector_seq_view.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/size_of.hpp>
#include <stan/math/prim/scal/meta/max_size.hpp>
#include <stan/math/prim/scal/meta/contains_nonconstant_struct.hpp>
#include <boost/random/normal_distribution.hpp>
#include <stan/math/prim/scal/fun/Phi.hpp>
#include <stan/math/prim/scal/prob/std_normal_lpdf.hpp>
#include <stan/math/prim/scal/fun/std_binormal_integral.hpp>
#include <stan/math/prim/scal/fun/inv.hpp>
#include <stan/math/prim/scal/fun/value_of_rec.hpp>
#include <stan/math/prim/scal/meta/max_size_mvt.hpp>
#include <stan/math/prim/mat/meta/value_type.hpp>

#include <cmath>
#include <limits>
#include <algorithm>

namespace stan {
namespace math {

/**
 * The log of the CDF of the standard bivariate normal (binormal) for the
 * specified vector(s) given the specified correlation(s). rho can
 * each be either a scalar or a vector. Any vector inputs must be the same
 * length.
 *
 * <p>The result log probability is defined to be the sum of the
 * log probabilities for each 2-vector observation/correlation
 * parameter pair.
 *
 * @param y (Sequence of) 2-vector(s).
 * @param rho (Sequence of) correlation parameter(s)
 * for the standard bivariate normal distribution.
 * @return The log of the product of the probabilities P(Y1 <= y[1], Y2 <= y[2]
 * | rho).
 * @throw std::domain_error if the rho is not between -1 and 1 or nan.
 * @tparam T_y Underlying type of random variable in sequence.
 * @tparam T_rho Underlying type of correlation in sequence.
 */

template <typename T_y, typename T_rho>
typename return_type<T_y, T_rho>::type std_binormal_lcdf(const T_y& y,
                                                         const T_rho& rho) {
  static const char* function = "std_binormal_lcdf";
  typedef
      typename stan::partials_return_type<T_y, T_rho>::type T_partials_return;
  typedef typename stan::math::value_type<T_y>::type T_y_child_type;

  using std::asin;
  using std::exp;
  using std::log;
  using std::max;
  using std::sqrt;
  using stan::math::value_of_rec;

  check_bounded(function, "Correlation parameter", rho, -1.0, 1.0);
  if (stan::is_vector_like<T_y_child_type>::value
      || stan::is_vector_like<T_rho>::value) {
    check_consistent_sizes(function, "Random variable", y,
                           "Correlation parameter", rho);
  }

  vector_seq_view<T_y> y_vec(y);
  const scalar_seq_view<T_rho> rho_vec(rho);
  size_t N_y = length_mvt(y);
  size_t N_rho = size_of(rho);
  size_t N = max(N_y, N_rho);

  operands_and_partials<T_y, T_rho> ops_partials(y, rho);

  T_partials_return cdf_log(0.0);
  if (unlikely(N == 0))
    return cdf_log;

  // check size consistency of all random variables y
  for (size_t i = 0; i < N_y; ++i) {
    check_size_match(function,
                     "Size of one of the vectors "
                     "of the location variable",
                     y_vec[i].size(),
                     "Size of another vector of "
                     "the location variable",
                     2);
  }

  for (size_t n = 0; n < N_y; ++n) {
    check_not_nan(function, "Random variable", y_vec[n]);
  }

  for (size_t n = 0; n < N; ++n) {
    const T_partials_return y1_dbl = value_of(y_vec[n][0]);
    const T_partials_return y2_dbl = value_of(y_vec[n][1]);
    const T_partials_return rho_dbl = value_of(rho_vec[n]);

    T_partials_return cdf_ = std_binormal_integral(y1_dbl, y2_dbl, rho_dbl);
    cdf_log += log(cdf_);

    if (contains_nonconstant_struct<T_y, T_rho>::value) {
      const T_partials_return inv_cdf_ = 1.0 / cdf_;
      const T_partials_return one_minus_rho_sq = (1 + rho_dbl) * (1 - rho_dbl);
      const T_partials_return sqrt_one_minus_rho_sq = sqrt(one_minus_rho_sq);
      const T_partials_return rho_times_y2 = rho_dbl * y2_dbl;
      const T_partials_return y1_minus_rho_times_y2 = y1_dbl - rho_times_y2;
      const bool y1_isfinite = std::isfinite(value_of_rec(y1_dbl));
      const bool y2_isfinite = std::isfinite(value_of_rec(y2_dbl));
      const bool y1_y2_arefinite = y1_isfinite && y2_isfinite;
      const bool rho_lt_one = fabs(rho_dbl) < 1;
      if (!is_constant_struct<T_y>::value) {
        ops_partials.edge1_.partials_vec_[n](0)
            += cdf_ > 0 && rho_lt_one && y1_y2_arefinite
               ? inv_cdf_ * exp(std_normal_lpdf(y1_dbl))
                * Phi((y2_dbl - rho_dbl * y1_dbl)
                      / sqrt_one_minus_rho_sq)
               : (!y2_isfinite && y1_isfinite && cdf_ > 0)
                 || (rho_dbl == 1 && y1_dbl < y2_dbl && cdf_ > 0)
                 || (rho_dbl == -1 && y2_dbl > -y1_dbl && cdf_ > 0)
               ? inv_cdf_ * exp(std_normal_lpdf(y1_dbl)) :
                cdf_ > 0 ? 0 : inv_cdf_;
        ops_partials.edge1_.partials_vec_[n](1)
            += cdf_ > 0 && rho_lt_one  && y1_y2_arefinite
               ? inv_cdf_ * exp(std_normal_lpdf(y2_dbl))
                         * Phi(y1_minus_rho_times_y2 / sqrt_one_minus_rho_sq)
               : (!y1_isfinite && y2_isfinite && cdf_ > 0)
                 || (rho_dbl == 1 && y2_dbl < y1_dbl && cdf_ > 0)
                 || (rho_dbl == -1 && y2_dbl > -y1_dbl && cdf_ > 0)
               ? inv_cdf_ * exp(std_normal_lpdf(y2_dbl)) :
                cdf_ > 0 ? 0 : inv_cdf_;
      }
      if (!is_constant_struct<T_rho>::value)
        ops_partials.edge2_.partials_[n]
            += cdf_ > 0 && y1_y2_arefinite && rho_lt_one
                   ? inv_cdf_ * 0.5 / (stan::math::pi() * sqrt_one_minus_rho_sq)
                         * exp(-0.5 / one_minus_rho_sq * y1_minus_rho_times_y2
                                   * y1_minus_rho_times_y2
                               - 0.5 * y2_dbl * y2_dbl) :
               cdf_ > 0 ? 0 : inv_cdf_;
    }
  }
  return ops_partials.build(cdf_log);
}

}  // namespace math
}  // namespace stan
#endif
