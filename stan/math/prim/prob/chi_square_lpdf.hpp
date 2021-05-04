#ifndef STAN_MATH_PRIM_PROB_CHI_SQUARE_LPDF_HPP
#define STAN_MATH_PRIM_PROB_CHI_SQUARE_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of a chi-squared density for y with the specified
 * degrees of freedom parameter.
 * The degrees of freedom parameter must be greater than 0.
 * y must be greater than or equal to 0.
 *
 \f{eqnarray*}{
 y &\sim& \chi^2_\nu \\
 \log (p (y \, |\, \nu)) &=& \log \left( \frac{2^{-\nu / 2}}{\Gamma (\nu / 2)}
 y^{\nu / 2 - 1} \exp^{- y / 2} \right) \\
 &=& - \frac{\nu}{2} \log(2) - \log (\Gamma (\nu / 2)) + (\frac{\nu}{2} - 1)
 \log(y) - \frac{y}{2} \\ & & \mathrm{ where } \; y \ge 0 \f}
 *
 * @tparam T_y type of scalar
 * @tparam T_dof type of degrees of freedom
 * @param y A scalar variable.
 * @param nu Degrees of freedom.
 * @throw std::domain_error if nu is not greater than or equal to 0
 * @throw std::domain_error if y is not greater than or equal to 0.
 */
template <bool propto, typename T_y, typename T_dof,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_dof>* = nullptr>
return_type_t<T_y, T_dof> chi_square_lpdf(const T_y& y, const T_dof& nu) {
  using T_partials_return = partials_return_t<T_y, T_dof>;
  using T_partials_array = Eigen::Array<T_partials_return, Eigen::Dynamic, 1>;
  using std::log;
  static const char* function = "chi_square_lpdf";
  using T_y_ref = ref_type_t<T_y>;
  using T_nu_ref = ref_type_t<T_dof>;
  check_consistent_sizes(function, "Random variable", y,
                         "Degrees of freedom parameter", nu);
  T_y_ref y_ref = y;
  T_nu_ref nu_ref = nu;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) nu_val = to_ref(as_value_column_array_or_scalar(nu_ref));

  check_nonnegative(function, "Random variable", y_val);
  check_positive_finite(function, "Degrees of freedom parameter", nu_val);

  if (size_zero(y, nu)) {
    return 0;
  }
  if (!include_summand<propto, T_y, T_dof>::value) {
    return 0;
  }

  size_t N = max_size(y, nu);
  const auto& log_y = to_ref_if<!is_constant_all<T_dof>::value>(log(y_val));
  const auto& half_nu = to_ref(0.5 * nu_val);

  T_partials_return logp(0);
  if (include_summand<propto, T_dof>::value) {
    logp -= sum(nu_val * HALF_LOG_TWO + lgamma(half_nu)) * N / size(nu);
  }
  logp += sum((half_nu - 1.0) * log_y);

  if (include_summand<propto, T_y>::value) {
    logp -= 0.5 * sum(y_val) * N / size(y);
  }

  operands_and_partials<T_y_ref, T_nu_ref> ops_partials(y_ref, nu_ref);
  if (!is_constant_all<T_y>::value) {
    ops_partials.edge1_.partials_ = (half_nu - 1.0) / y_val - 0.5;
  }
  if (!is_constant_all<T_dof>::value) {
    if (is_vector<T_dof>::value) {
      ops_partials.edge2_.partials_ = forward_as<T_partials_array>(
          (log_y - digamma(half_nu)) * 0.5 - HALF_LOG_TWO);
    } else {
      ops_partials.edge2_.partials_[0]
          = sum(log_y - digamma(half_nu)) * 0.5 - HALF_LOG_TWO * N;
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_dof>
inline return_type_t<T_y, T_dof> chi_square_lpdf(const T_y& y,
                                                 const T_dof& nu) {
  return chi_square_lpdf<false>(y, nu);
}

}  // namespace math
}  // namespace stan
#endif
