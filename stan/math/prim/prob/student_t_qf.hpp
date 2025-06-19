#ifndef STAN_MATH_PRIM_PROB_STUDENT_T_QF_HPP
#define STAN_MATH_PRIM_PROB_STUDENT_T_QF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/inv_inc_beta.hpp>

namespace stan {
namespace math {

double student_t_qf(const double p, const double nu, const double mu,
                    const double sigma) {
  static constexpr const char* function = "student_t_qf";
  check_nonnegative(function, "Degrees of freedom parameter", nu);
  check_positive(function, "Scale parameter", sigma);
  check_bounded(function, "Probability parameter", p, 0.0, 1.0);

  if (p == 0.0) {
    return NEGATIVE_INFTY;
  } else if (p == 1.0) {
    return INFTY;
  } else if (p == 0.5) {
    return mu;
  }

  const double p_val_flip = p < 0.5 ? p : 1.0 - p;
  const double p_sign = p < 0.5 ? -1.0 : 1.0;
  const auto& ibeta_arg = inv_inc_beta(0.5 * nu, 0.5, 2 * p_val_flip);

  return mu + p_sign * sigma * math::sqrt(nu)
          * math::sqrt(-1.0 + 1.0 / ibeta_arg);
}
}
}

#endif
