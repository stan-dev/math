#ifndef STAN_MATH_PRIM_MAT_PROB_BINORMAL_COPULA_CDF_HPP
#define STAN_MATH_PRIM_MAT_PROB_BINORMAL_COPULA_CDF_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/mat/meta/operands_and_partials.hpp>
#include <stan/math/prim/mat/meta/vector_seq_view.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
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
#include <stan/math/prim/scal/fun/inv_Phi.hpp>
#include <stan/math/prim/scal/fun/std_binormal_integral.hpp>

#include <cmath>
#include <limits>
#include <algorithm>

namespace stan {
namespace math {

/**
 * The CDF of the bivariate normal (binormal) copula for the
 * specified vector given the specified correlation.
 *
 * Reference for the gradients of the function:
 * Meyer, Christian. "The bivariate normal copula." Communications
 * in Statistics-Theory and Methods 42, no. 13 (2013): 2402-2422.
 * https://arxiv.org/pdf/0912.2816.pdf
 *
 * @param u 2-vector.
 * @param rho correlation parameter
 * for the standard bivariate normal distribution.
 * @return The probability P(U1 <= u[1], U2 <= u[2] | rho).
 * @throw std::domain_error if the rho is not between -1 and 1 or nan.
 * @throw std::domain_error if elements of u aren't between 0 and 1 or nan.
 * @throw std::invalid_argument u isn't size 2.
 * @tparam T_u Type of random variable.
 * @tparam T_rho Type of correlation.
 */

template <typename T_u, typename T_rho>
typename return_type<T_u, T_rho>::type binormal_copula_cdf(const T_u& u,
                                                           const T_rho& rho) {
  static const char* function = "binormal_copula_cdf";
  typedef
      typename stan::partials_return_type<T_u, T_rho>::type T_partials_return;

  using std::asin;
  using std::exp;
  using std::log;
  using std::max;
  using std::sqrt;

  check_bounded(function, "Correlation parameter", rho, -1.0, 1.0);
  check_bounded(function, "Random variable", u, 0.0, 1.0);
  check_size_match(function,
                   "Size of the vector "
                   "of the random variable",
                   u.size(),
                   "the dimension of the "
                   "the support of the distribution",
                   2);

  operands_and_partials<T_u, T_rho> ops_partials(u, rho);

  const T_partials_return u1_dbl = value_of(u[0]);
  const T_partials_return u2_dbl = value_of(u[1]);
  const T_partials_return rho_dbl = value_of(rho);
  const T_partials_return z1_dbl = inv_Phi(u1_dbl);
  const T_partials_return z2_dbl = inv_Phi(u2_dbl);

  T_partials_return cdf_ = std_binormal_integral(z1_dbl, z2_dbl, rho_dbl);

  if (contains_nonconstant_struct<T_u, T_rho>::value) {
    const T_partials_return one_minus_rho_sq = (1 + rho_dbl) * (1 - rho_dbl);
    const T_partials_return sqrt_one_minus_rho_sq = sqrt(one_minus_rho_sq);
    const T_partials_return rho_times_z2 = rho_dbl * z2_dbl;
    const T_partials_return z1_minus_rho_times_z2 = z1_dbl - rho_times_z2;
    if (!is_constant_struct<T_u>::value) {
      ops_partials.edge1_.partials_vec_[0](0)
          += Phi((z2_dbl - rho_dbl * z1_dbl) / sqrt_one_minus_rho_sq);
      ops_partials.edge1_.partials_vec_[0](1)
          += Phi(z1_minus_rho_times_z2 / sqrt_one_minus_rho_sq);
    }
    if (!is_constant_struct<T_rho>::value)
      ops_partials.edge2_.partials_[0]
          += 0.5 / (stan::math::pi() * sqrt_one_minus_rho_sq)
             * exp(-0.5 / one_minus_rho_sq * z1_minus_rho_times_z2
                       * z1_minus_rho_times_z2
                   - 0.5 * z2_dbl * z2_dbl);
  }
  return ops_partials.build(cdf_);
}

}  // namespace math
}  // namespace stan
#endif
