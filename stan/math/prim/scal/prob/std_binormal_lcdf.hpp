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
#include <stan/math/prim/scal/fun/Phi.hpp>
#include <stan/math/prim/scal/prob/std_normal_lpdf.hpp>
#include <stan/math/prim/scal/fun/binormal_integral_owens.hpp>
#include <stan/math/prim/scal/fun/inv.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

template <typename T_y_1, typename T_y_2, typename T_rho>
typename return_type<T_y_1, T_y_2, T_rho>::type std_binormal_lcdf(
    const T_y_1& y_1, const T_y_2& y_2, const T_rho& rho) {
  static const char* function = "std_binormal_lcdf";
  typedef typename stan::partials_return_type<T_y_1, T_y_2, T_rho>::type
      T_partials_return;

  using std::exp;
  using std::sqrt;
  using std::log;
  using std::asin;

  T_partials_return cdf_log(0.0);
  if (size_zero(y_1, y_2, rho))
    return cdf_log;

  check_not_nan(function, "Random variable 1", y_1);
  check_not_nan(function, "Random variable 2", y_2);
  check_bounded(function, "Correlation parameter", rho, -1.0, 1.0);
  check_consistent_sizes(function, "Random variable 1", y_1, "Random variable 2",
                         y_2, "Correlation parameter", rho);

  operands_and_partials<T_y_1, T_y_2, T_rho> ops_partials(y_1, y_2, rho);

  scalar_seq_view<T_y_1> y_1_vec(y_1);
  scalar_seq_view<T_y_2> y_2_vec(y_2);
  scalar_seq_view<T_rho> rho_vec(rho);
  size_t N = max_size(y_1, y_2, rho);

  for (size_t n = 0; n < N; ++n) {
    const T_partials_return y_1_dbl = value_of(y_1_vec[n]);
    const T_partials_return y_2_dbl = value_of(y_2_vec[n]);
    const T_partials_return rho_dbl = value_of(rho_vec[n]);

    T_partials_return cdf_ = binormal_integral_owens(y_1_dbl, y_2_dbl, rho_dbl);
      cdf_log += log(cdf_);

    if (contains_nonconstant_struct<T_y_1, T_y_2, T_rho>::value) {
      const T_partials_return inv_cdf_ = cdf_ > 0 ? inv(cdf_) : 
        std::numeric_limits<double>::infinity();
      const T_partials_return one_minus_rho_sq = (1 + rho_dbl) * (1 - rho_dbl);
      const T_partials_return sqrt_one_minus_rho_sq = sqrt(one_minus_rho_sq);
      const T_partials_return rho_times_y_2 = rho_dbl * y_2_dbl;
      const T_partials_return y_1_minus_rho_times_y_2 = y_1_dbl - rho_times_y_2;
      if (!is_constant_struct<T_y_1>::value)
        ops_partials.edge1_.partials_[n] 
          += cdf_ > 0 && cdf_ < 1 ? inv_cdf_ * exp(std_normal_lpdf(y_1_dbl))
              * Phi((y_2_dbl - rho_dbl * y_1_dbl) / sqrt_one_minus_rho_sq)
              : cdf_ > 0 ? 1 : 0;
      if (!is_constant_struct<T_y_2>::value)
        ops_partials.edge2_.partials_[n] 
          += cdf_ > 0 && cdf_ < 1 ? inv_cdf_ * exp(std_normal_lpdf(y_2_dbl))
              * Phi(y_1_minus_rho_times_y_2 / sqrt_one_minus_rho_sq)
              : cdf_ > 0 ? 1 : 0;
      if (!is_constant_struct<T_rho>::value)
        ops_partials.edge3_.partials_[n]
            += cdf_ > 0 && cdf_ < 1 ? inv_cdf_ * 0.5 / (stan::math::pi() * sqrt_one_minus_rho_sq)
            * exp(-0.5 / one_minus_rho_sq * y_1_minus_rho_times_y_2 *  y_1_minus_rho_times_y_2
                  -0.5 * y_2_dbl * y_2_dbl)
              : cdf_ > 0 ? 1 : 0;
    }
  }
  return ops_partials.build(cdf_log);
}

}  // namespace math
}  // namespace stan
#endif
