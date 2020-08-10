#ifndef STAN_MATH_PRIM_PROB_STUDENT_T_CDF_HPP
#define STAN_MATH_PRIM_PROB_STUDENT_T_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/beta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/grad_reg_inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_dof, typename T_loc, typename T_scale>
return_type_t<T_y, T_dof, T_loc, T_scale> student_t_cdf(const T_y& y,
                                                        const T_dof& nu,
                                                        const T_loc& mu,
                                                        const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_dof, T_loc, T_scale>;
  using std::exp;
  using std::pow;
  static const char* function = "student_t_cdf";
  check_not_nan(function, "Random variable", y);
  check_positive_finite(function, "Degrees of freedom parameter", nu);
  check_finite(function, "Location parameter", mu);
  check_positive_finite(function, "Scale parameter", sigma);

  if (size_zero(y, nu, mu, sigma)) {
    return 1.0;
  }

  T_partials_return P(1.0);
  operands_and_partials<T_y, T_dof, T_loc, T_scale> ops_partials(y, nu, mu,
                                                                 sigma);
  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_dof> nu_vec(nu);
  scalar_seq_view<T_loc> mu_vec(mu);
  scalar_seq_view<T_scale> sigma_vec(sigma);
  size_t N = max_size(y, nu, mu, sigma);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < stan::math::size(y); i++) {
    if (value_of(y_vec[i]) == NEGATIVE_INFTY) {
      return ops_partials.build(0.0);
    }
  }

  T_partials_return digammaHalf = 0;

  VectorBuilder<!is_constant_all<T_dof>::value, T_partials_return, T_dof>
      digamma_vec(size(nu));
  VectorBuilder<!is_constant_all<T_dof>::value, T_partials_return, T_dof>
      digammaNu_vec(size(nu));
  VectorBuilder<!is_constant_all<T_dof>::value, T_partials_return, T_dof>
      digammaNuPlusHalf_vec(size(nu));

  if (!is_constant_all<T_dof>::value) {
    digammaHalf = digamma(0.5);

    for (size_t i = 0; i < stan::math::size(nu); i++) {
      const T_partials_return nu_dbl = value_of(nu_vec[i]);

      digammaNu_vec[i] = digamma(0.5 * nu_dbl);
      digammaNuPlusHalf_vec[i] = digamma(0.5 + 0.5 * nu_dbl);
    }
  }

  for (size_t n = 0; n < N; n++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (value_of(y_vec[n]) == INFTY) {
      continue;
    }

    const T_partials_return sigma_inv = 1.0 / value_of(sigma_vec[n]);
    const T_partials_return t
        = (value_of(y_vec[n]) - value_of(mu_vec[n])) * sigma_inv;
    const T_partials_return nu_dbl = value_of(nu_vec[n]);
    const T_partials_return q = nu_dbl / (t * t);
    const T_partials_return r = 1.0 / (1.0 + q);
    const T_partials_return J = 2 * r * r * q / t;
    const T_partials_return betaNuHalf = beta(0.5, 0.5 * nu_dbl);
    double zJacobian = t > 0 ? -0.5 : 0.5;

    if (q < 2) {
      T_partials_return z
          = inc_beta(0.5 * nu_dbl, (T_partials_return)0.5, 1.0 - r);
      const T_partials_return Pn = t > 0 ? 1.0 - 0.5 * z : 0.5 * z;
      const T_partials_return d_ibeta
          = pow(r, -0.5) * pow(1.0 - r, 0.5 * nu_dbl - 1) / betaNuHalf;

      P *= Pn;

      if (!is_constant_all<T_y>::value) {
        ops_partials.edge1_.partials_[n]
            += -zJacobian * d_ibeta * J * sigma_inv / Pn;
      }
      if (!is_constant_all<T_dof>::value) {
        T_partials_return g1 = 0;
        T_partials_return g2 = 0;

        grad_reg_inc_beta(g1, g2, 0.5 * nu_dbl, (T_partials_return)0.5, 1.0 - r,
                          digammaNu_vec[n], digammaHalf,
                          digammaNuPlusHalf_vec[n], betaNuHalf);

        ops_partials.edge2_.partials_[n]
            += zJacobian * (d_ibeta * (r / t) * (r / t) + 0.5 * g1) / Pn;
      }

      if (!is_constant_all<T_loc>::value) {
        ops_partials.edge3_.partials_[n]
            += zJacobian * d_ibeta * J * sigma_inv / Pn;
      }
      if (!is_constant_all<T_scale>::value) {
        ops_partials.edge4_.partials_[n]
            += zJacobian * d_ibeta * J * sigma_inv * t / Pn;
      }

    } else {
      T_partials_return z
          = 1.0 - inc_beta((T_partials_return)0.5, 0.5 * nu_dbl, r);

      zJacobian *= -1;

      const T_partials_return Pn = t > 0 ? 1.0 - 0.5 * z : 0.5 * z;

      T_partials_return d_ibeta
          = pow(1.0 - r, 0.5 * nu_dbl - 1) * pow(r, -0.5) / betaNuHalf;

      P *= Pn;

      if (!is_constant_all<T_y>::value) {
        ops_partials.edge1_.partials_[n]
            += zJacobian * d_ibeta * J * sigma_inv / Pn;
      }
      if (!is_constant_all<T_dof>::value) {
        T_partials_return g1 = 0;
        T_partials_return g2 = 0;

        grad_reg_inc_beta(g1, g2, (T_partials_return)0.5, 0.5 * nu_dbl, r,
                          digammaHalf, digammaNu_vec[n],
                          digammaNuPlusHalf_vec[n], betaNuHalf);

        ops_partials.edge2_.partials_[n]
            += zJacobian * (-d_ibeta * (r / t) * (r / t) + 0.5 * g2) / Pn;
      }
      if (!is_constant_all<T_loc>::value) {
        ops_partials.edge3_.partials_[n]
            += -zJacobian * d_ibeta * J * sigma_inv / Pn;
      }
      if (!is_constant_all<T_scale>::value) {
        ops_partials.edge4_.partials_[n]
            += -zJacobian * d_ibeta * J * sigma_inv * t / Pn;
      }
    }
  }

  if (!is_constant_all<T_y>::value) {
    for (size_t n = 0; n < stan::math::size(y); ++n) {
      ops_partials.edge1_.partials_[n] *= P;
    }
  }
  if (!is_constant_all<T_dof>::value) {
    for (size_t n = 0; n < stan::math::size(nu); ++n) {
      ops_partials.edge2_.partials_[n] *= P;
    }
  }
  if (!is_constant_all<T_loc>::value) {
    for (size_t n = 0; n < stan::math::size(mu); ++n) {
      ops_partials.edge3_.partials_[n] *= P;
    }
  }
  if (!is_constant_all<T_scale>::value) {
    for (size_t n = 0; n < stan::math::size(sigma); ++n) {
      ops_partials.edge4_.partials_[n] *= P;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
