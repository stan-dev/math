#ifndef STAN_MATH_PRIM_PROB_LOGLOGISTIC_CDF_HPP
#define STAN_MATH_PRIM_PROB_LOGLOGISTIC_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>
#include <iostream>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of the loglogistic density for the specified scalar(s)
 * given the specified scales(s) and shape(s). y, alpha, or beta
 * can each be either a scalar or a vector. Any vector inputs
 * must be the same length.
 *
 * <p>The result log probability is defined to be the sum of the
 * log probabilities for each observation/scale/shape triple.
 *
 * @tparam T_y type of scalar.
 * @tparam T_scale type of scale parameter.
 * @tparam T_shape type of shape parameter.
 * @param y (Sequence of) scalar(s).
 * @param alpha (Sequence of) scale parameter(s)
 * for the loglogistic distribution.
 * @param beta (Sequence of) shape parameter(s) for the
 * loglogistic distribution.
 * @return The log of the product of the densities.
 * @throw std::domain_error if any of the inputs are not positive and finite.
 */
template <typename T_y, typename T_scale, typename T_shape,
           require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
           T_y, T_scale, T_shape>* = nullptr>
return_type_t<T_y, T_scale, T_shape> loglogistic_cdf(const T_y& y,
                                                     const T_scale& alpha,
                                                     const T_shape& beta) {
  using T_partials_return = partials_return_t<T_y, T_scale, T_shape>;
  using T_y_ref = ref_type_t<T_y>;
  using T_alpha_ref = ref_type_t<T_scale>;
  using T_beta_ref = ref_type_t<T_shape>;
  using std::pow;
  static const char* function = "loglogistic_cdf";
  check_consistent_sizes(function, "Random variable", y, "Scale parameter",
                         alpha, "Shape parameter", beta);
  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) alpha_val = to_ref(as_value_column_array_or_scalar(alpha_ref));
  decltype(auto) beta_val = to_ref(as_value_column_array_or_scalar(beta_ref));


  check_nonnegative(function, "Random variable", y_val);
  check_positive_finite(function, "Scale parameter", alpha_val);
  check_positive_finite(function, "Shape parameter", beta_val);

  if (size_zero(y, alpha, beta)) {
    return 1.0;
  }

  operands_and_partials<T_y_ref, T_alpha_ref, T_beta_ref> ops_partials(
      y_ref, alpha_ref, beta_ref);

  const auto& cdf_elementwise = to_ref_if<!is_constant_all<T_y, T_scale,
    T_shape>::value>(1 / (1 + pow(alpha_val / y_val, beta_val)));

  T_partials_return cdf(1.0);
  cdf *= prod(cdf_elementwise);

  if (!is_constant_all<T_y>::value) {
    const auto& y_deriv = (pow(alpha_val, beta_val) * beta_val *
      pow(y_val, -beta_val - 1)) / pow(1 + pow(alpha_val / y_val, beta_val), 2);
    ops_partials.edge1_.partials_ = (y_deriv / cdf_elementwise) * cdf;
    // std::cout << "Partial 1: " << y_deriv << std::endl;
  }
  if (!is_constant_all<T_scale>::value) {
    const auto& alpha_deriv = (-pow(y_val, -beta_val) * beta_val *
      pow(alpha_val, beta_val - 1)) / pow(1 + pow(alpha_val / y_val, beta_val), 2);
    ops_partials.edge2_.partials_ = (alpha_deriv / cdf_elementwise) * cdf;
    // std::cout << "Partial 2: " << alpha_deriv << std::endl;
    // std::cout << "Partial 2: " << (alpha_deriv == 0) << std::endl;
  }
  if (!is_constant_all<T_shape>::value) {
    const auto& beta_deriv = (-pow(alpha_val / y_val, beta_val) *
      log(alpha_val / y_val)) /
      pow(1 + pow(alpha_val / y_val, beta_val), 2);
    ops_partials.edge3_.partials_ = (beta_deriv / cdf_elementwise) * cdf;
    // Is it the product here, rather than sum, if it is a scalar???


    // std::cout << "Partial 3: " << beta_deriv << std::endl;
  }
  if (!is_constant_all<T_shape>::value) {
    for (size_t n = 0; n < stan::math::size(beta); ++n) {
      // std::cout << "ii1" << std::endl;
      if (ops_partials.edge3_.partials_[n] != ops_partials.edge3_.partials_[n]) {
        ops_partials.edge3_.partials_[n] = 0;
        // std::cout << "ii2" << std::endl;
      }
      // ops_partials.edge3_.partials_[n] *= cdf;
    }
    // std::cout << "Partial 3: " << ops_partials.edge3_.partials_ << std::endl;
  }
  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
