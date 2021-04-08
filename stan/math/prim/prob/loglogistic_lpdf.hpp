#ifndef STAN_MATH_PRIM_PROB_LOGLOGISTIC_LPDF_HPP
#define STAN_MATH_PRIM_PROB_LOGLOGISTIC_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

#include <iostream>
#include <typeinfo>

namespace stan {
namespace math {

  // Kako testiramo? Ali se testi nekje drugje pokličejo?
  // Izgleda, kot da se kliče tisti class nekje.


// Logistic(y|mu, sigma)    [sigma > 0]
template <bool propto, typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
return_type_t<T_y, T_loc, T_scale> loglogistic_lpdf(const T_y& y, const T_loc& mu,
                                                 const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_mu_ref = ref_type_if_t<!is_constant<T_loc>::value, T_loc>;
  using T_sigma_ref = ref_type_if_t<!is_constant<T_scale>::value, T_scale>;
  using std::pow;
  static const char* function = "loglogistic_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);

  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));
  // Eigen::Matrix<double, 1, 4> tmppow;
  // tmppow = pow(y_val / mu_val, sigma_val);
  // T_partials_return logp = sum(tmppow);
  // std::cout << std::endl;
  // std::cout << typeid(y_val).name() << std::endl;
  // std::cout << std::endl;
  // std::cout << type_name<decltype(y_val)>() << '\n';
  // std::cout << std::endl;


  check_finite(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive_finite(function, "Scale parameter", sigma_val);

  if (size_zero(y, mu, sigma)) {
    return 0.0;
  }
  if (!include_summand<propto, T_y, T_loc, T_scale>::value) {
    return 0.0;
  }

  // Here we create some class that will be put onto the stack I guess.
  operands_and_partials<T_y_ref, T_mu_ref, T_sigma_ref> ops_partials(
      y_ref, mu_ref, sigma_ref);

  const auto& inv_sigma
      = to_ref_if<!is_constant_all<T_y, T_loc, T_scale>::value>(inv(sigma_val));
  const auto& y_minus_mu
      = to_ref_if<!is_constant_all<T_scale>::value>(y_val - mu_val);
  const auto& y_minus_mu_div_sigma = to_ref(y_minus_mu * inv_sigma);

  const auto& y_div_mu = to_ref(y_val / mu_val);

  size_t N = max_size(y, mu, sigma);
  // The propto part I guess.
  //T_partials_return logp = (log(sigma_val) - log(mu_val) + (sigma_val - 1) * (log(y_val) - log(mu_val))) - 2 * log1p(std::pow((y_val / mu_val), sigma_val));
  // T_partials_return logp = - 2 * std::log1p(std::pow((y_val / mu_val), sigma_val));
  // T_partials_return logp = std::pow((y_val / mu_val), sigma_val);
  // T_partials_return logp = std::pow((1 / 2), 3);
  // T_partials_return logp = std::pow((y_val / (double)mu_val), sigma_val); // Something like this. Is it OK?
  // T_partials_return logp = sum(std::pow((y_val / mu_val), sigma_val));
  // T_partials_return logp = sum(pow((y_val / mu_val), sigma_val));
  // T_partials_return logp = sum(exp((y_val / mu_val)));
  // T_partials_return logp = sum(pow(y_div_mu, sigma_val));
  // T_partials_return logp = pow(((double)y_val / (double)mu_val), (double)sigma_val);
  // T_partials_return logp = 2;
  // T_partials_return logp = sum(pow((y_val / mu_val), sigma_val));
  // T_partials_return logp = sum(pow(y_div_mu, sigma_val));

  //T_partials_return logp = pow((1 / 2), 3);
  T_partials_return logp = sum((log(sigma_val) - log(mu_val) + (sigma_val - 1) *
    (log(y_val) - log(mu_val))) -
    2 * log1p(pow((y_val / mu_val), sigma_val)));

  // OK this works. I guess I can add the derivatives and later see,
  // which parts would make sense to save as separate computations.
  // Also TODO: dividing integers!




  // real lalpha = log(alpha);
// real numerator = log(beta) - lalpha + (beta - 1) * (log(y) - lalpha);
// return numerator - 2 * (log1p( (y/alpha)^beta ));
  // T_partials_return logp = -sum(y_minus_mu_div_sigma)
  //                          - 2.0 * sum(log1p(exp(-y_minus_mu_div_sigma)));
  // if (include_summand<propto, T_scale>::value) { // Do we include normalization? Is this it?
  //   logp -= sum(log(sigma_val)) * N / size(sigma);
  // }

  // Here come the derivatives I think. We divide everything up by
  // being constant or not.
  if (!is_constant_all<T_y, T_scale>::value) {
    const auto& exp_y_minus_mu_div_sigma = exp(y_minus_mu_div_sigma);
    const auto& y_deriv = to_ref_if<(!is_constant_all<T_scale>::value
                                     && !is_constant_all<T_y>::value)>(
        (2 / (1 + exp_y_minus_mu_div_sigma) - 1) * inv_sigma);
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_ = y_deriv;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_ = (-y_deriv * y_minus_mu - 1) * inv_sigma;
    }
  }
  if (!is_constant_all<T_loc>::value) {
    const auto& exp_mu_div_sigma = to_ref(exp(mu_val * inv_sigma));
    ops_partials.edge2_.partials_
        = (1
           - 2 * exp_mu_div_sigma / (exp_mu_div_sigma + exp(y_val * inv_sigma)))
          * inv_sigma;
  }
  return ops_partials.build(logp);
}

// Ahaa, this is just the default call of loglogistic with the
// normalization constant.
template <typename T_y, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_loc, T_scale> loglogistic_lpdf(const T_y& y,
                                                        const T_loc& mu,
                                                        const T_scale& sigma) {
  return loglogistic_lpdf<false>(y, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
