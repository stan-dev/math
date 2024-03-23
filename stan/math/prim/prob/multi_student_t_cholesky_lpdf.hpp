#ifndef STAN_MATH_PRIM_PROB_MULTI_STUDENT_T_CHOLESKY_LPDF_HPP
#define STAN_MATH_PRIM_PROB_MULTI_STUDENT_T_CHOLESKY_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/columns_dot_self.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/dot_product.hpp>
#include <stan/math/prim/fun/dot_self.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/is_inf.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/max_size_mvt.hpp>
#include <stan/math/prim/fun/mdivide_left_tri.hpp>
#include <stan/math/prim/fun/mdivide_right_tri.hpp>
#include <stan/math/prim/fun/size_mvt.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/vector_seq_view.hpp>
#include <stan/math/prim/prob/multi_normal_cholesky_lpdf.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>
#include <cstdlib>
#include <iostream>

namespace stan {
namespace math {
/** \ingroup multivar_dists
 * The log of the multivariate student t density for the given y, mu,
 * nu, and a Cholesky factor L of the scale matrix.
 * Sigma = LL', a square, semi-positive definite matrix.
 *
 * This version of the function is vectorized on y and mu.
 *
 * @param y scalar vector of random variates
 * @param nu scalar degrees of freedom
 * @param mu location vector
 * @param L cholesky decomposition of the scale matrix
 * of the multivariate student t distribution.
 * @return The log of the multivariate student t density.
 * @throw std::domain_error if LL' is not square, not symmetric,
 * or not semi-positive definite.
 * @tparam T_y Type of scalar.
 * @tparam T_dof Type of scalar.
 * @tparam T_loc Type of location.
 * @tparam T_covar Type of scale.
 * @return log probability of the multivariate student t distribution.
 */
template <
    bool propto, typename T_y, typename T_dof, typename T_loc, typename T_covar,
    require_any_not_vector_vt<is_stan_scalar, T_y, T_dof, T_loc>* = nullptr,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
        T_y, T_loc, T_covar>* = nullptr>
return_type_t<T_y, T_dof, T_loc, T_covar> multi_student_t_cholesky_lpdf(
    const T_y& y, const T_dof& nu, const T_loc& mu, const T_covar& L) {
  static constexpr const char* function = "multi_student_t_cholesky";
  using T_covar_elem = typename scalar_type<T_covar>::type;
  using T_return = return_type_t<T_y, T_dof, T_loc, T_covar>;
  using T_partials_return = partials_return_t<T_y, T_dof, T_loc, T_covar>;
  using matrix_partials_t
      = Eigen::Matrix<T_partials_return, Eigen::Dynamic, Eigen::Dynamic>;
  using vector_partials_t = Eigen::Matrix<T_partials_return, Eigen::Dynamic, 1>;
  using row_vector_partials_t
      = Eigen::Matrix<T_partials_return, 1, Eigen::Dynamic>;
  using T_y_ref = ref_type_t<T_y>;
  using T_nu_ref = ref_type_t<T_dof>;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_L_ref = ref_type_t<T_covar>;

  check_consistent_sizes_mvt(function, "y", y, "mu", mu);
  size_t num_y = size_mvt(y);
  size_t num_mu = size_mvt(mu);

  if (num_y == 0 || num_mu == 0) {
    return 0;
  }

  T_y_ref y_ref = y;
  T_nu_ref nu_ref = nu;
  T_mu_ref mu_ref = mu;
  T_L_ref L_ref = L;
  vector_seq_view<T_y_ref> y_vec(y_ref);
  vector_seq_view<T_mu_ref> mu_vec(mu_ref);
  const size_t size_vec = max_size_mvt(y_ref, mu_ref);

  check_not_nan(function, "Degrees of freedom parameter", nu_ref);
  check_positive(function, "Degrees of freedom parameter", nu_ref);
  check_finite(function, "Degrees of freedom parameter", nu_ref);
  check_cholesky_factor(function, "scale parameter", L_ref);

  const int size_y = y_vec[0].size();
  const int size_mu = mu_vec[0].size();
  const int num_dims = L.rows();

  if (unlikely(num_dims == 0)) {
    return T_return(0);
  }

  for (size_t i = 1, size_mvt_y = num_y; i < size_mvt_y; i++) {
    check_size_match(
        function, "Size of one of the vectors of the random variable",
        y_vec[i].size(), "Size of another vector of the random variable",
        y_vec[i - 1].size());
  }

  for (size_t i = 1, size_mvt_mu = num_mu; i < size_mvt_mu; i++) {
    check_size_match(function,
                     "Size of one of the vectors "
                     "of the location variable",
                     mu_vec[i].size(),
                     "Size of another vector of "
                     "the location variable",
                     mu_vec[i - 1].size());
  }

  check_size_match(function, "Size of random variable", size_mu,
                   "rows of scale parameter", L.rows());
  check_size_match(function, "Size of random variable", size_y,
                   "size of location parameter", size_mu);
  check_size_match(function, "Size of random variable", size_y,
                   "rows of scale parameter", L.rows());
  check_size_match(function, "Size of random variable", size_y,
                   "columns of scale parameter", L.cols());

  for (size_t i = 0; i < size_vec; i++) {
    check_finite(function, "Location parameter", mu_vec[i]);
    check_not_nan(function, "Random variable", y_vec[i]);
  }

  auto ops_partials = make_partials_propagator(y_ref, nu_ref, mu_ref, L_ref);

  T_partials_return lp(0);

  if (include_summand<propto>::value) {
    lp += -0.5 * num_dims * LOG_PI * size_vec;
  }

  if (include_summand<propto, T_y, T_dof, T_loc, T_covar_elem>::value) {
    T_partials_return nu_val = value_of(nu_ref);
    T_partials_return inv_nu = inv(nu_val);
    T_partials_return nu_plus_dims = nu_val + num_dims;
    matrix_partials_t L_val = value_of(L_ref);
    matrix_partials_t L_deriv;
    const auto& half_nu
        = to_ref_if<include_summand<propto, T_dof>::value>(0.5 * nu_val);
    const auto& digamma_vals = to_ref_if<!is_constant<T_dof>::value>(
        digamma(half_nu + 0.5 * num_dims) - digamma(half_nu));

    if (include_summand<propto, T_dof>::value) {
      lp += lgamma(0.5 * nu_plus_dims) * size_vec;
      lp += -lgamma(0.5 * nu_val) * size_vec;
      lp += -(0.5 * num_dims) * log(nu_val) * size_vec;
    }

    if (include_summand<propto, T_covar_elem>::value) {
      lp += -sum(log(L_val.diagonal())) * size_vec;
    }

    T_partials_return sum_lp_vec(0.0);
    row_vector_partials_t half(size_y);
    vector_partials_t y_val_minus_mu_val(size_y);
    vector_partials_t scaled_diff(size_y);

    for (size_t i = 0; i < size_vec; i++) {
      decltype(auto) y_val = as_value_column_vector_or_scalar(y_vec[i]);
      decltype(auto) mu_val = as_value_column_vector_or_scalar(mu_vec[i]);
      y_val_minus_mu_val = eval(y_val - mu_val);

      half = mdivide_left_tri<Eigen::Lower>(L_val, y_val_minus_mu_val)
                 .transpose();

      scaled_diff = mdivide_right_tri<Eigen::Lower>(half, L_val).transpose();

      T_partials_return dot_half = dot_self(half);

      if (!is_constant_all<T_dof>::value) {
        T_partials_return G = dot_product(scaled_diff, y_val_minus_mu_val);
        partials<1>(ops_partials)[i] += 0.5
                                        * (digamma_vals - log1p(G * inv_nu)
                                           + (G - num_dims) / (G + nu_val));
      }

      scaled_diff *= nu_plus_dims / (dot_half + nu_val);

      if (!is_constant_all<T_y>::value) {
        partials_vec<0>(ops_partials)[i] += -scaled_diff;
      }

      if (!is_constant_all<T_loc>::value) {
        partials_vec<2>(ops_partials)[i] += scaled_diff;
      }

      if (!is_constant_all<T_covar_elem>::value) {
        if (i == 0) {
          L_deriv
              = (scaled_diff * half).template triangularView<Eigen::Lower>();
        } else {
          L_deriv
              += (scaled_diff * half).template triangularView<Eigen::Lower>();
        }
      }

      sum_lp_vec += log1p(dot_half * inv_nu);
    }

    if (!is_constant_all<T_covar_elem>::value) {
      L_deriv.diagonal().array() -= size_vec / L_val.diagonal().array();
      partials<3>(ops_partials) += L_deriv;
    }
    lp += -0.5 * nu_plus_dims * sum_lp_vec;
  }
  return ops_partials.build(lp);
}

/** \ingroup multivar_dists
 * The log of the multivariate student t density for the given y, mu,
 * nu, and a Cholesky factor L of the scale matrix.
 * Sigma = LL', a square, semi-positive definite matrix.
 *
 * This version of the function is vectorized on y and mu.
 *
 * @param y A scalar vector.
 * @param nu A scalar degrees of freedom.
 * @param mu A location vector.
 * @param L The Cholesky decomposition of the scale matrix
 * of the multivariate student t distribution.
 * @return The log of the multivariate student t density.
 * @throw std::domain_error if LL' is not square, not symmetric,
 * or not semi-positive definite.
 * @tparam T_y Type of scalar.
 * @tparam T_dof Type of scalar.
 * @tparam T_loc Type of location.
 * @tparam T_covar Type of scale.
 * @return log probability of the multivariate student t distribution.
 */
template <bool propto, typename T_y, typename T_dof, typename T_loc,
          typename T_covar,
          require_all_vector_vt<is_stan_scalar, T_y, T_dof, T_loc>* = nullptr,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_dof, T_loc, T_covar>* = nullptr>
return_type_t<T_y, T_dof, T_loc, T_covar> multi_student_t_cholesky_lpdf(
    const T_y& y, const T_dof& nu, const T_loc& mu, const T_covar& L) {
  static const char* function = "multi_student_t_cholesky";
  using T_covar_elem = typename scalar_type<T_covar>::type;
  using Eigen::Matrix;
  using T_return = return_type_t<T_y, T_dof, T_loc, T_covar>;
  using T_partials_return = partials_return_t<T_y, T_dof, T_loc, T_covar>;
  using matrix_partials_t
      = Eigen::Matrix<T_partials_return, Eigen::Dynamic, Eigen::Dynamic>;
  using vector_partials_t = Eigen::Matrix<T_partials_return, Eigen::Dynamic, 1>;
  using row_vector_partials_t
      = Eigen::Matrix<T_partials_return, 1, Eigen::Dynamic>;
  using T_y_ref = ref_type_t<T_y>;
  using T_nu_ref = ref_type_t<T_dof>;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_L_ref = ref_type_t<T_covar>;

  T_y_ref y_ref = y;
  T_nu_ref nu_ref = nu;
  T_mu_ref mu_ref = mu;
  T_L_ref L_ref = L;
  vector_partials_t y_val = as_value_column_vector_or_scalar(y_ref);
  vector_partials_t mu_val = as_value_column_vector_or_scalar(mu_ref);

  const int size_y = y_ref.size();
  const int size_mu = mu_ref.size();

  if (unlikely(size_y == 0)) {
    return T_return(0);
  }

  check_not_nan(function, "Degrees of freedom parameter", nu_ref);
  check_positive(function, "Degrees of freedom parameter", nu_ref);
  check_finite(function, "Degrees of freedom parameter", nu_ref);

  check_size_match(function, "Size of random variable", size_mu,
                   "rows of scale parameter", L.rows());
  check_size_match(function, "Size of random variable", size_y,
                   "size of location parameter", size_mu);
  check_size_match(function, "Size of random variable", size_y,
                   "rows of scale parameter", L.rows());
  check_size_match(function, "Size of random variable", size_y,
                   "columns of scale parameter", L.cols());

  check_finite(function, "Location parameter", mu_val);
  check_not_nan(function, "Random variable", y_val);

  check_cholesky_factor(function, "scale parameter", L_ref);

  T_partials_return lp(0);

  auto ops_partials = make_partials_propagator(y_ref, nu_ref, mu_ref, L_ref);

  if (include_summand<propto>::value) {
    lp += -0.5 * size_y * LOG_PI;
  }

  if (include_summand<propto, T_y, T_dof, T_loc, T_covar_elem>::value) {
    T_partials_return nu_val = value_of(nu_ref);
    T_partials_return inv_nu = inv(nu_val);
    T_partials_return nu_plus_dims = nu_val + size_y;
    vector_partials_t y_val_minus_mu_val = eval(y_val - mu_val);

    matrix_partials_t L_val = value_of(L_ref);
    row_vector_partials_t half
        = mdivide_left_tri<Eigen::Lower>(L_val, y_val_minus_mu_val).transpose();
    vector_partials_t scaled_diff
        = mdivide_right_tri<Eigen::Lower>(half, L_val).transpose();
    T_partials_return dot_half = dot_self(half);

    if (!is_constant_all<T_dof>::value) {
      T_partials_return half_nu = 0.5 * nu_val;
      T_partials_return digamma_vals
          = digamma(half_nu + 0.5 * size_y) - digamma(half_nu);
      T_partials_return G = dot_product(scaled_diff, y_val_minus_mu_val);

      partials<1>(ops_partials)
          += 0.5
             * (digamma_vals - log1p(G * inv_nu) + (G - size_y) / (G + nu_val));
    }

    if (include_summand<propto, T_dof>::value) {
      lp += lgamma(0.5 * (nu_val + size_y));
      lp += -lgamma(0.5 * nu_val);
      lp += -0.5 * size_y * log(nu_val);
    }

    if (include_summand<propto, T_covar_elem>::value) {
      lp += -sum(log(L_val.diagonal()));
    }

    if (!is_constant_all<T_y, T_loc, T_covar_elem>::value) {
      T_partials_return scale_val = nu_plus_dims / (dot_half + nu_val);

      if (!is_constant_all<T_y>::value) {
        partials<0>(ops_partials) += -scaled_diff * scale_val;
      }
      if (!is_constant_all<T_loc>::value) {
        partials<2>(ops_partials) += scaled_diff * scale_val;
      }
      if (!is_constant_all<T_covar_elem>::value) {
        matrix_partials_t L_deriv
            = (scaled_diff * half).template triangularView<Eigen::Lower>();
        L_deriv.diagonal().array() -= 1 / L_val.diagonal().array();
        edge<3>(ops_partials).partials_ += L_deriv;
      }
    }

    lp += -0.5 * nu_plus_dims * log1p(dot_half * inv_nu);
  }

  return ops_partials.build(lp);
}

template <typename T_y, typename T_dof, typename T_loc, typename T_covar>
inline return_type_t<T_y, T_dof, T_loc, T_covar> multi_student_t_cholesky_lpdf(
    const T_y& y, const T_dof& nu, const T_loc& mu, const T_covar& L) {
  return multi_student_t_cholesky_lpdf<false>(y, nu, mu, L);
}

}  // namespace math
}  // namespace stan
#endif
