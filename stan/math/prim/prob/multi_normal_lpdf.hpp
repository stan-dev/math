#ifndef STAN_MATH_PRIM_PROB_MULTI_NORMAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_MULTI_NORMAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/cholesky_decompose.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/max_size_mvt.hpp>
#include <stan/math/prim/fun/mdivide_left_tri.hpp>
#include <stan/math/prim/fun/mdivide_right_tri.hpp>
#include <stan/math/prim/prob/multi_normal_cholesky_lpdf.hpp>
#include <stan/math/prim/fun/trace_inv_quad_form_ldlt.hpp>
#include <stan/math/prim/fun/vector_seq_view.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>

namespace stan {
namespace math {

template <bool propto, typename T_y, typename T_loc, typename T_covar,
          require_any_not_vector_vt<is_stan_scalar, T_y, T_loc>* = nullptr,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_covar>* = nullptr>
return_type_t<T_y, T_loc, T_covar> multi_normal_lpdf(const T_y& y,
                                                     const T_loc& mu,
                                                     const T_covar& Sigma) {
  static const char* function = "multi_normal_lpdf";
  using T_y_ref = ref_type_t<T_y>;
  using T_mu_ref = ref_type_t<T_loc>;

  check_consistent_sizes_mvt(function, "y", y, "mu", mu);
  size_t number_of_y = size_mvt(y);
  size_t number_of_mu = size_mvt(mu);
  if (number_of_y == 0 || number_of_mu == 0) {
    return 0;
  }

  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;

  vector_seq_view<T_y_ref> y_vec(y_ref);
  vector_seq_view<T_mu_ref> mu_vec(mu_ref);
  const size_t size_vec = max_size_mvt(y, mu);

  const int size_y = y_vec[0].size();
  const int size_mu = mu_vec[0].size();

  // check size consistency of all random variables y
  for (size_t i = 1, size_mvt_y = size_mvt(y); i < size_mvt_y; i++) {
    check_size_match(function,
                     "Size of one of the vectors of "
                     "the random variable",
                     y_vec[i].size(),
                     "Size of the first vector of the "
                     "random variable",
                     size_y);
  }
  // check size consistency of all means mu
  for (size_t i = 1, size_mvt_mu = size_mvt(mu); i < size_mvt_mu; i++) {
    check_size_match(function,
                     "Size of one of the vectors of "
                     "the location variable",
                     mu_vec[i].size(),
                     "Size of the first vector of the "
                     "location variable",
                     size_mu);
  }

  check_size_match(function, "Size of random variable", size_y,
                   "size of location parameter", size_mu);
  check_size_match(function, "Size of random variable", size_y,
                   "rows of covariance parameter", Sigma.rows());
  check_size_match(function, "Size of random variable", size_y,
                   "columns of covariance parameter", Sigma.cols());

  for (size_t i = 0; i < size_vec; i++) {
    check_finite(function, "Location parameter", mu_vec[i]);
    check_not_nan(function, "Random variable", y_vec[i]);
  }

  if (unlikely(size_y == 0)) {
    return T_return(0);
  }

  auto L = cholesky_decompose(Sigma);
  return multi_normal_cholesky_lpdf<propto>(y, mu, L);
}

template <bool propto, typename T_y, typename T_loc, typename T_covar,
          require_all_vector_vt<is_stan_scalar, T_y, T_loc>* = nullptr,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_covar>* = nullptr>
return_type_t<T_y, T_loc, T_covar> multi_normal_lpdf(const T_y& y,
                                                     const T_loc& mu,
                                                     const T_covar& Sigma) {
  static const char* function = "multi_normal_lpdf";
  using T_y_ref = ref_type_t<T_y>;
  using T_mu_ref = ref_type_t<T_loc>;

  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;

  decltype(auto) y_val = as_value_column_vector_or_scalar(y_ref);
  decltype(auto) mu_val = as_value_column_vector_or_scalar(mu_ref);

  const int size_y = y_ref.size();
  const int size_mu = mu_ref.size();

  check_size_match(function, "Size of random variable", size_y,
                   "size of location parameter", size_mu);
  check_size_match(function, "Size of random variable", size_y,
                   "rows of covariance parameter", Sigma.rows());
  check_size_match(function, "Size of random variable", size_y,
                   "columns of covariance parameter", Sigma.cols());

  check_finite(function, "Location parameter", mu_val);
  check_not_nan(function, "Random variable", y_val);

  if (unlikely(size_y == 0)) {
    return T_return(0);
  }

  auto L = cholesky_decompose(Sigma);
  return multi_normal_cholesky_lpdf<propto>(y, mu, L);
}

template <typename T_y, typename T_loc, typename T_covar>
inline return_type_t<T_y, T_loc, T_covar> multi_normal_lpdf(const T_y& y,
                                                            const T_loc& mu,
                                                            const T_covar& L) {
  return multi_normal_lpdf<false>(y, mu, L);
}

}  // namespace math
}  // namespace stan
#endif
