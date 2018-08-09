#ifndef STAN_MATH_PRIM_SCAL_PROB_BETA_PROPORTION_LCCDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BETA_PROPORTION_LCCDF_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_less_or_equal.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/digamma.hpp>
#include <stan/math/prim/scal/fun/lbeta.hpp>
#include <stan/math/prim/scal/meta/contains_nonconstant_struct.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/scal/fun/grad_reg_inc_beta.hpp>
#include <stan/math/prim/scal/fun/inc_beta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the beta log complementary cumulative distribution function
 * for the given probability and parameterized with location and
 * precision parameters. Given matching containers returns the log sum
 * of probabilities.
 *
 * @tparam T_y type of probability parameter
 * @tparam T_loc type of location parameter
 * @tparam T_prec type of precision parameter
 * @param y probability parameter
 * @param p location parameter
 * @param c precision parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if p is outside (0, 1)
 * @throw std::domain_error if c is nonpositive
 * @throw std::domain_error if y is not a valid probability
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_y, typename T_loc, typename T_prec>
typename return_type<T_y, T_loc, T_prec>::type beta_proportion_lccdf(
    const T_y& y, const T_loc& p, const T_prec& c) {

  typedef
    typename stan::partials_return_type<T_y, T_loc, T_prec>::type
    T_partials_return;

  static const char* function = "beta_proportion_lccdf";

  if (size_zero(y, p, c))
    return 0.0;

  using boost::math::tools::promote_args;

  T_partials_return ccdf_log(0.0);

  check_positive(function, "Location parameter", p);
  check_less_or_equal(function, "Location parameter", p, 1.0);
  check_positive_finite(function, "Precision parameter", c);
  check_not_nan(function, "Random variable", y);
  check_nonnegative(function, "Random variable", y);
  check_less_or_equal(function, "Random variable", y, 1.0);
  check_consistent_sizes(function, "Random variable", y,
                         "Location parameter", p,
                         "Precision parameter", c);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_loc> p_vec(p);
  scalar_seq_view<T_prec> c_vec(c);
  size_t N = max_size(y, p, c);

  operands_and_partials<T_y, T_loc, T_prec> ops_partials(y, p, c);

  using std::exp;
  using std::exp;
  using std::log;
  using std::pow;

  VectorBuilder<contains_nonconstant_struct<T_loc, T_prec>::value,
                T_partials_return, T_loc, T_prec>
      digamma_pc(max_size(p, c));
  VectorBuilder<contains_nonconstant_struct<T_loc, T_prec>::value,
                T_partials_return, T_loc, T_prec>
      digamma_c_1m_p(max_size(p, c));
  VectorBuilder<contains_nonconstant_struct<T_loc, T_prec>::value,
                T_partials_return, T_loc, T_prec>
    digamma_c(max_size(p, c));

  if (contains_nonconstant_struct<T_loc, T_prec>::value) {
    for (size_t i = 0; i < N; i++) {
      const T_partials_return pc_dbl = value_of(p_vec[i]) * value_of(c_vec[i]);
      const T_partials_return c_1m_p_dbl
        = value_of(c_vec[i]) * (1.0 - value_of(p_vec[i]));

      digamma_pc[i] = digamma(pc_dbl);
      digamma_c_1m_p[i] = digamma(c_1m_p_dbl);
      digamma_c[i] = digamma(value_of(c_vec[i]));
    }
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return p_dbl = value_of(p_vec[n]);
    const T_partials_return c_dbl = value_of(c_vec[n]);
    const T_partials_return pc_dbl = p_dbl * c_dbl;
    const T_partials_return c_1m_p_dbl = c_dbl * (1.0 - p_dbl);
    const T_partials_return betafunc_dbl = exp(lbeta(pc_dbl, c_1m_p_dbl));
    const T_partials_return Pn = 1.0 - inc_beta(pc_dbl, c_1m_p_dbl, y_dbl);

    ccdf_log += log(Pn);

    if (!is_constant_struct<T_y>::value)
      ops_partials.edge1_.partials_[n] -= pow(1 - y_dbl, c_1m_p_dbl - 1)
                                          * pow(y_dbl, pc_dbl - 1)
                                          / betafunc_dbl / Pn;

    T_partials_return g1 = 0;
    T_partials_return g2 = 0;

    if (contains_nonconstant_struct<T_loc, T_prec>::value) {
      grad_reg_inc_beta(g1, g2, pc_dbl, c_1m_p_dbl, y_dbl,
                        digamma_pc[n], digamma_c_1m_p[n],
                        digamma_c[n], betafunc_dbl);
    }
    if (!is_constant_struct<T_loc>::value)
      ops_partials.edge2_.partials_[n] -= c_dbl * (g1 - g2) / Pn;
    if (!is_constant_struct<T_prec>::value)
      ops_partials.edge3_.partials_[n]
        -= (g1 * p_dbl + g2 * (1.0 - p_dbl)) / Pn;
  }
  return ops_partials.build(ccdf_log);
}

}  // namespace math
}  // namespace stan
#endif
