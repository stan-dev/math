#ifndef STAN_MATH_PRIM_SCAL_PROB_STD_BINORMAL_LCDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_STD_BINORMAL_LCDF_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/max_size.hpp>
#include <stan/math/prim/scal/meta/contains_nonconstant_struct.hpp>
#include <boost/random/normal_distribution.hpp>
#include <stan/math/prim/scal/fun/inv.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

/**
 * The log of the pdf of the standard bivariate normal (binormal) for the
 * specified scalar(s) given the specified correlation(s). y1, y2, or rho can
 * each be either a scalar or a vector. Any vector inputs must be the same
 * length.
 *
 * <p>The result log density is defined to be the sum of the
 * log densities for each observation 1/observation 2/correlation 
 * parameter triple.
 *
 * @param y1 (Sequence of) scalar(s).
 * @param y2 (Sequence of) scalar(s).
 * @param rho (Sequence of) correlation parameter(s)
 * for the standard bivariate normal distribution.
 * @return The log of the product of the probabilities P(Y1 <= y1, Y2 <= y2 | rho).
 * @throw std::domain_error if the rho is not between -1 and 1.
 * @tparam T_y1 Underlying type of random variable 1 in sequence.
 * @tparam T_y2 Underlying type of random variable 2 in sequence.
 * @tparam T_rho Underlying type of correlation in sequence.
 */

template <bool propto, typename T_y1, typename T_y2, typename T_rho>
typename return_type<T_y1, T_y2, T_rho>::type std_binormal_lpdf(
    const T_y1& y1, const T_y2& y2, const T_rho& rho) {
  static const char* function = "std_binormal_lpdf";
  typedef typename stan::partials_return_type<T_y1, T_y2, T_rho>::type
      T_partials_return;

  using std::log;

  T_partials_return logp(0.0);
  if (size_zero(y1, y2, rho))
    return logp;

  check_not_nan(function, "Random variable 1", y1);
  check_not_nan(function, "Random variable 2", y2);
  check_bounded(function, "Correlation parameter", rho, -1.0, 1.0);
  check_consistent_sizes(function, "Random variable 1", y1, "Random variable 2",
                         y2, "Correlation parameter", rho);

  if (!include_summand<propto, T_y1, T_y2, T_rho>::value)
    return 0.0;

  operands_and_partials<T_y1, T_y2, T_rho> ops_partials(y1, y2, rho);

  scalar_seq_view<T_y1> y1_vec(y1);
  scalar_seq_view<T_y2> y2_vec(y2);
  scalar_seq_view<T_rho> rho_vec(rho);
  size_t N = max_size(y1, y2, rho);

  VectorBuilder<true, T_partials_return, T_rho> one_minus_rho_sq_(length(rho));
  VectorBuilder<include_summand<propto, T_rho>::value, T_partials_return,
                T_rho>
      log_one_minus_rho_sq(length(sigma));

  for (size_t i = 0; i < length(sigma); ++i) {
    T_partials_return rho_dbl = value_of(rho_vec[i]);
    one_minus_rho_sq[i] = (1 + rho_dbl) * (1 - rho_dbl);
    if (include_summand<propto, T_rho>::value)
      log_one_minus_rho_sq[i] = log(one_minus_rho_sq[i]);
  }

  for (size_t n = 0; n < N; ++n) {
    const T_partials_return y1_dbl = value_of(y1_vec[n]);
    const T_partials_return y2_dbl = value_of(y2_vec[n]);
    const T_partials_return rho_dbl = value_of(rho_vec[n]);
    const T_partials_return y1_sq = y1_dbl * y1_dbl;
    const T_partials_return y2_sq = y2_dbl * y2_dbl;
    const T_partials_return two_rho_y1_y2 = 2 * rho_dbl * y1_dbl * y2_dbl;
    const T_partials_return quad_form = y1_sq + y2_sq - two_rho_y1_y2;
    const T_partials_return quad_form_over_one_minus_rho_sq = quad_form / one_minus_rho_sq[n];

    static double NEGATIVE_HALF = -0.5;

    if (include_summand<propto, T_rho>::value)
      logp += NEGATIVE_HALF * log_one_minus_rho_sq[n];

    if (include_summand<propto>::value)
      logp += NEG_LOG_SQRT_TWO_PI;

    if (include_summand<propto, T_y1, T_y2>::value)
      logp += NEGATIVE_HALF * quad_form_over_one_minus_rho_sq;

    if (!is_constant_struct<T_y1>::value)
      ops_partials.edge1_.partials_[n] 
        += (-y1_dbl + y2_dbl * rho_dbl) / one_minus_rho_sq[n];
    if (!is_constant_struct<T_y2>::value)
      ops_partials.edge2_.partials_[n] 
        += (-y2_dbl + y1_dbl * rho_dbl) / one_minus_rho_sq[n];
    if (!is_constant_struct<T_rho>::value)
      ops_partials.edge3_.partials_[n]
          += rho_dbl / one_minus_rho_sq[n] * (1.0 - quad_form_over_one_minus_rho_sq);
  }
  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
