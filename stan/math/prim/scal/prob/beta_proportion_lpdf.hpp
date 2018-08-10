#ifndef STAN_MATH_PRIM_SCAL_PROB_BETA_PROPORTION_LPDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BETA_PROPORTION_LPDF_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_less_or_equal.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/log1m.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/digamma.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <stan/math/prim/scal/meta/contains_nonconstant_struct.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * The log of the beta density parameterized for the location and
 * precision. y, p (location), or c (precision) can each either be
 * scalar or a vector.  Any vector inputs must be the same length.
 *
 * <p> The result log probability is defined to be the sum of
 * the log probabilities for each observation/p/c triple.
 *
 * Prior location, p, must be contained in (0, 1).  Prior precision
 * must be positive.
 *
 * @param y (Sequence of) scalar(s).
 * @param p (Sequence of) prior location(s).
 * @param c (Sequence of) prior precision(s).
 * @return The log of the product of densities.
 * @tparam T_y Type of scalar outcome.
 * @tparam T_loc Type of prior location.
 * @tparam T_prec Type of prior precision.
 */
template <bool propto, typename T_y, typename T_loc, typename T_prec>
typename return_type<T_y, T_loc, T_prec>::type beta_proportion_lpdf(
    const T_y& y, const T_loc& p, const T_prec& c) {
  static const char* function = "beta_proportion_lpdf";

  typedef typename stan::partials_return_type<T_y, T_loc, T_prec>::type
      T_partials_return;

  using stan::is_constant_struct;
  using std::log;

  if (size_zero(y, p, c))
    return 0.0;

  T_partials_return logp(0.0);

  check_positive(function, "Location parameter", p);
  check_less_or_equal(function, "Location parameter", p, 1.0);
  check_positive_finite(function, "Precision parameter", c);
  check_not_nan(function, "Random variable", y);
  check_nonnegative(function, "Random variable", y);
  check_less_or_equal(function, "Random variable", y, 1.0);
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         p, "Precision parameter", c);

  if (!include_summand<propto, T_y, T_loc, T_prec>::value)
    return 0.0;

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_loc> p_vec(p);
  scalar_seq_view<T_prec> c_vec(c);
  size_t N = max_size(y, p, c);

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    if (y_dbl < 0 || y_dbl > 1)
      return LOG_ZERO;
  }

  operands_and_partials<T_y, T_loc, T_prec> ops_partials(y, p, c);

  VectorBuilder<include_summand<propto, T_y, T_loc, T_prec>::value,
                T_partials_return, T_y>
      log_y(length(y));
  VectorBuilder<include_summand<propto, T_y, T_loc, T_prec>::value,
                T_partials_return, T_y>
      log1m_y(length(y));

  for (size_t n = 0; n < length(y); n++) {
    if (include_summand<propto, T_y, T_loc, T_prec>::value) {
      log_y[n] = log(value_of(y_vec[n]));
      log1m_y[n] = log1m(value_of(y_vec[n]));
    }
  }

  VectorBuilder<include_summand<propto, T_loc, T_prec>::value,
                T_partials_return, T_loc, T_prec>
      lgamma_pc(max_size(p, c));
  VectorBuilder<contains_nonconstant_struct<T_loc, T_prec>::value,
                T_partials_return, T_loc, T_prec>
      digamma_pc(max_size(p, c));

  for (size_t n = 0; n < max_size(p, c); n++) {
    const T_partials_return pc = value_of(p_vec[n]) * value_of(c_vec[n]);
    if (include_summand<propto, T_loc, T_prec>::value)
      lgamma_pc[n] = lgamma(pc);
    if (contains_nonconstant_struct<T_loc, T_prec>::value)
      digamma_pc[n] = digamma(pc);
  }

  VectorBuilder<include_summand<propto, T_loc, T_prec>::value,
                T_partials_return, T_loc, T_prec>
      lgamma_c_1m_p(max_size(p, c));
  VectorBuilder<contains_nonconstant_struct<T_loc, T_prec>::value,
                T_partials_return, T_loc, T_prec>
      digamma_c_1m_p(max_size(p, c));

  for (size_t n = 0; n < max_size(p, c); n++) {
    const T_partials_return c_1m_p
        = value_of(c_vec[n]) * (1.0 - value_of(p_vec[n]));
    if (include_summand<propto, T_loc, T_prec>::value)
      lgamma_c_1m_p[n] = lgamma(c_1m_p);
    if (contains_nonconstant_struct<T_loc, T_prec>::value)
      digamma_c_1m_p[n] = digamma(c_1m_p);
  }

  VectorBuilder<include_summand<propto, T_prec>::value, T_partials_return,
                T_prec>
      lgamma_c(length(c));
  VectorBuilder<!is_constant_struct<T_prec>::value, T_partials_return, T_prec>
      digamma_c(length(c));

  for (size_t n = 0; n < length(c); n++) {
    if (include_summand<propto, T_prec>::value)
      lgamma_c[n] = lgamma(value_of(c_vec[n]));
    if (!is_constant_struct<T_prec>::value)
      digamma_c[n] = digamma(value_of(c_vec[n]));
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return p_dbl = value_of(p_vec[n]);
    const T_partials_return c_dbl = value_of(c_vec[n]);

    if (include_summand<propto, T_prec>::value)
      logp += lgamma_c[n];
    if (include_summand<propto, T_loc, T_prec>::value)
      logp -= lgamma_pc[n] + lgamma_c_1m_p[n];
    if (include_summand<propto, T_y, T_loc, T_prec>::value) {
      const T_partials_return pc_dbl = p_dbl * c_dbl;
      logp += (pc_dbl - 1.0) * log_y[n] + (c_dbl - pc_dbl - 1.0) * log1m_y[n];
    }

    if (!is_constant_struct<T_y>::value) {
      const T_partials_return pc_dbl = p_dbl * c_dbl;
      ops_partials.edge1_.partials_[n]
          += (pc_dbl - 1.0) / y_dbl + (c_dbl - pc_dbl - 1.0) / (y_dbl - 1.0);
    }
    if (!is_constant_struct<T_loc>::value)
      ops_partials.edge2_.partials_[n]
          += c_dbl
             * (digamma_c_1m_p[n] - digamma_pc[n] + log_y[n] - log1m_y[n]);
    if (!is_constant_struct<T_prec>::value)
      ops_partials.edge3_.partials_[n]
          += digamma_c[n] + p_dbl * (log_y[n] - digamma_pc[n])
             + (1.0 - p_dbl) * (log1m_y[n] - digamma_c_1m_p[n]);
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_prec>
inline typename return_type<T_y, T_loc, T_prec>::type beta_proportion_lpdf(
    const T_y& y, const T_loc& p, const T_prec& c) {
  return beta_proportion_lpdf<false>(y, p, c);
}

}  // namespace math
}  // namespace stan
#endif
