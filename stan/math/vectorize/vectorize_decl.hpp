
#ifndef STAN_MATH_VECTORIZE_VECTORIZE_DECL_HPP
#define STAN_MATH_VECTORIZE_VECTORIZE_DECL_HPP
#include <stan/math/prim/meta.hpp>
namespace stan {
namespace math {

template <typename Container, require_not_stan_scalar_t<Container> *  = nullptr, require_not_container_st<std::is_arithmetic, Container> *  = nullptr, require_not_var_matrix_t<Container> *  = nullptr, require_not_nonscalar_prim_or_rev_kernel_expression_t<Container> *  = nullptr>
 inline  auto square(const Container & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr>
 inline  auto square(const Container & x);

template <typename T, require_not_container_st<std::is_arithmetic, T> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr>
 inline  auto inv(const T & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<Container> *  = nullptr>
 inline  auto inv(const Container & x);

template <typename Container, require_not_container_st<std::is_arithmetic, Container> *  = nullptr, require_not_var_matrix_t<Container> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<Container> *  = nullptr, require_not_stan_scalar_t<Container> *  = nullptr>
 inline  auto fabs(const Container & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr>
 inline  auto fabs(const Container & x);

template <typename Container, require_not_container_st<std::is_arithmetic, Container> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<Container> *  = nullptr, require_not_var_matrix_t<Container> *  = nullptr>
 inline  auto floor(const Container & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr, require_not_var_matrix_t<Container> *  = nullptr>
 inline  auto floor(const Container & x);

template <typename Container, require_not_container_st<std::is_arithmetic, Container> *  = nullptr, require_not_var_matrix_t<Container> *  = nullptr, require_not_stan_scalar_t<Container> *  = nullptr>
 inline  auto abs(const Container & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr>
 inline  auto abs(const Container & x);

template <typename Container, require_not_container_st<std::is_arithmetic, Container> *  = nullptr, require_not_var_matrix_t<Container> *  = nullptr, require_not_nonscalar_prim_or_rev_kernel_expression_t<Container> *  = nullptr>
 inline  auto log(const Container & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr>
 inline  auto log(const Container & x);

template <typename T, require_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr, require_not_var_matrix_t<T> *  = nullptr>
 inline  auto log1p(const T & x);

template <typename T, require_not_var_matrix_t<T> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr>
 inline  auto log1m(const T & x);

template <typename Container, require_not_container_st<std::is_arithmetic, Container> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<Container> *  = nullptr, require_not_var_matrix_t<Container> *  = nullptr>
 inline  auto sqrt(const Container & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr, require_not_var_matrix_t<Container> *  = nullptr>
 inline  auto sqrt(const Container & x);

template <typename Container, require_not_container_st<std::is_arithmetic, Container> *  = nullptr, require_not_nonscalar_prim_or_rev_kernel_expression_t<Container> *  = nullptr, require_not_var_matrix_t<Container> *  = nullptr>
 inline  auto exp(const Container & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr>
 inline  auto exp(const Container & x);

template <typename Container, require_not_container_st<std::is_arithmetic, Container> *  = nullptr, require_not_var_matrix_t<Container> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<Container> *  = nullptr>
 inline  auto tanh(const Container & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr>
 inline  auto tanh(const Container & x);

template <bool Jacobian, typename T, require_std_vector_t<T> *  = nullptr>
 inline  auto positive_constrain(const T & x, return_type_t<T> & lp);

template <typename T, require_not_var_matrix_t<T> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr>
 inline  auto inv_logit(const T & x);

template <typename T, require_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr, require_not_var_matrix_t<T> *  = nullptr>
 inline  auto log1p_exp(const T & x);

template <typename Container, require_not_container_st<std::is_arithmetic, Container> *  = nullptr, require_not_var_matrix_t<Container> *  = nullptr, require_not_stan_scalar_t<Container> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<Container> *  = nullptr>
 inline  auto inv_sqrt(const Container & x);

template <typename Container, require_not_var_matrix_t<Container> *  = nullptr, require_container_st<std::is_arithmetic, Container> *  = nullptr>
 inline  auto inv_sqrt(const Container & x);

template <bool Jacobian, typename T, require_std_vector_t<T> *  = nullptr>
 inline  auto sum_to_zero_constrain(const T & y, return_type_t<T> & lp);

template <typename Container, require_not_container_st<std::is_arithmetic, Container> *  = nullptr, require_not_var_matrix_t<Container> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<Container> *  = nullptr>
 inline  auto cosh(const Container & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr>
 inline  auto cosh(const Container & x);

template <typename Container, require_not_container_st<std::is_arithmetic, Container> *  = nullptr, require_not_var_matrix_t<Container> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<Container> *  = nullptr>
 inline  auto cos(const Container & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr>
 inline  auto cos(const Container & x);

template <typename T, require_not_var_matrix_t<T> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr>
 inline  auto asinh(const T & x);

template <typename Container, require_not_container_st<std::is_arithmetic, Container> *  = nullptr, require_not_var_matrix_t<Container> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<Container> *  = nullptr>
 inline  auto asin(const Container & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr>
 inline  auto asin(const Container & x);

template <typename Container, require_not_container_st<std::is_arithmetic, Container> *  = nullptr, require_not_var_matrix_t<Container> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<Container> *  = nullptr>
 inline  auto acos(const Container & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr>
 inline  auto acos(const Container & x);

template <typename T, require_not_var_matrix_t<T> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr>
 inline  auto acosh(const T & x);

template <typename T, require_not_var_matrix_t<T> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr>
 inline  auto atanh(const T & x);

template <typename Container, require_not_container_st<std::is_arithmetic, Container> *  = nullptr, require_not_var_matrix_t<Container> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<Container> *  = nullptr>
 inline  auto atan(const Container & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr>
 inline  auto atan(const Container & x);

template <typename T, require_not_var_matrix_t<T> *  = nullptr, require_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr>
 inline  auto lgamma(const T & x);

template <typename T, require_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr, require_not_var_matrix_t<T> *  = nullptr>
 inline  auto digamma(const T & x);

template <typename T, require_not_var_matrix_t<T> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr>
 inline  auto cbrt(const T & x);

template <typename Container, require_not_container_st<std::is_arithmetic, Container> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<Container> *  = nullptr>
 inline  auto ceil(const Container & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr, require_not_var_matrix_t<Container> *  = nullptr>
 inline  auto ceil(const Container & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr>
 inline  auto inv_square(const Container & x);

template <typename T, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr, require_not_var_matrix_t<T> *  = nullptr>
 inline  auto erf(const T & x);

template <typename T, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr, require_not_var_matrix_t<T> *  = nullptr>
 inline  auto erfc(const T & x);

template <typename T, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr, require_not_var_matrix_t<T> *  = nullptr>
 inline  auto exp2(const T & x);

template <typename T, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr, require_not_var_matrix_t<T> *  = nullptr>
 inline  auto expm1(const T & x);

template <typename Container, require_not_container_st<std::is_arithmetic, Container> *  = nullptr, require_not_var_matrix_t<Container> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<Container> *  = nullptr>
 inline  auto sinh(const Container & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr>
 inline  auto sinh(const Container & x);

template <typename T, require_not_container_st<std::is_arithmetic, T> *  = nullptr, require_not_var_matrix_t<T> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr>
 inline  auto sin(const T & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr>
 inline  auto sin(const Container & x);

template <typename T, require_container_t<T> *  = nullptr>
 inline  auto sign(const T & x);

template <typename T, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr, require_not_var_matrix_t<T> *  = nullptr>
 inline  auto tgamma(const T & x);

template <typename T, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr, require_not_var_matrix_t<T> *  = nullptr>
 inline  auto Phi(const T & x);

template <typename T, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr, require_not_var_matrix_t<T> *  = nullptr>
 inline  auto inv_Phi(const T & x);

template <typename Container, require_not_container_st<std::is_arithmetic, Container> *  = nullptr, require_not_var_matrix_t<Container> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<Container> *  = nullptr>
 inline  auto inv_cloglog(const Container & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr>
 inline  auto inv_cloglog(const Container & x);

template <typename T, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr, require_not_var_matrix_t<T> *  = nullptr, require_not_arithmetic_t<T> *  = nullptr>
 inline  auto inv_erfc(const T & x);

template <typename T, require_not_stan_scalar_t<T> *  = nullptr, require_not_var_matrix_t<T> *  = nullptr>
 inline  auto lambert_w0(const T & x);

template <typename T, require_not_stan_scalar_t<T> *  = nullptr, require_not_var_matrix_t<T> *  = nullptr>
 inline  auto lambert_wm1(const T & x);

template <typename Container, require_not_var_matrix_t<Container> *  = nullptr, require_not_container_st<std::is_arithmetic, Container> *  = nullptr, require_not_nonscalar_prim_or_rev_kernel_expression_t<Container> *  = nullptr>
 inline  auto log10(const Container & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr>
 inline  auto log10(const Container & x);

template <typename T, require_not_var_matrix_t<T> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr>
 inline  auto log1m_exp(const T & x);

template <typename T, require_not_var_matrix_t<T> *  = nullptr, require_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr>
 inline  typename apply_scalar_unary<log1m_inv_logit_fun, T>::return_t log1m_inv_logit(const T & x);

template <typename T, require_not_var_matrix_t<T> *  = nullptr, require_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr>
 inline  auto log2(const T & x);

template <typename T, require_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr, require_not_var_matrix_t<T> *  = nullptr>
 inline  auto log_inv_logit(const T & x);

template <typename Container, require_st_arithmetic<Container> *  = nullptr, require_container_t<Container> *  = nullptr>
 inline  auto log_softmax(const Container & x);

template <typename Container, require_not_container_st<std::is_arithmetic, Container> *  = nullptr, require_not_var_matrix_t<Container> *  = nullptr, require_not_nonscalar_prim_or_rev_kernel_expression_t<Container> *  = nullptr>
 inline  auto logit(const Container & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr>
 inline  auto logit(const Container & x);

template <typename T>
 inline  auto minus(const std::vector<T> & x);

template <typename T, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr, require_not_var_matrix_t<T> *  = nullptr>
 inline  auto Phi_approx(const T & x);

template <typename Container, require_not_container_st<std::is_arithmetic, Container> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<Container> *  = nullptr>
 inline  auto round(const Container & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr>
 inline  auto round(const Container & x);

template <typename T, require_container_t<T> *  = nullptr>
 inline  auto step(const T & x);

template <typename Container, require_not_container_st<std::is_arithmetic, Container> *  = nullptr, require_not_var_matrix_t<Container> *  = nullptr, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<Container> *  = nullptr>
 inline  auto tan(const Container & x);

template <typename Container, require_container_st<std::is_arithmetic, Container> *  = nullptr>
 inline  auto tan(const Container & x);

template <typename Container, require_std_vector_st<std::is_arithmetic, Container> *  = nullptr>
 inline  auto to_int(const Container & x);

template <typename T, require_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr>
 inline  auto trigamma(const T & x);

template <typename T, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr>
 inline  auto trunc(const T & x);

template <bool Jacobian, typename T, require_std_vector_t<T> *  = nullptr>
 inline  auto cholesky_corr_constrain(const T & y, int K, return_type_t<T> & lp);

template <typename T, require_std_vector_t<T> *  = nullptr>
  auto cholesky_corr_free(const T & x);

template <bool Jacobian, typename T, require_std_vector_t<T> *  = nullptr>
 inline  auto cholesky_factor_constrain(const T & x, int M, int N, return_type_t<T> & lp);

template <typename T, require_std_vector_t<T> *  = nullptr>
  auto cholesky_factor_free(const T & x);

template <bool Jacobian, typename T, require_std_vector_t<T> *  = nullptr>
 inline  auto corr_matrix_constrain(const T & y, int K, return_type_t<T> & lp);

template <typename T, require_std_vector_t<T> *  = nullptr>
  auto corr_matrix_free(const T & x);

template <bool Jacobian, typename T, require_std_vector_t<T> *  = nullptr>
 inline  auto cov_matrix_constrain(const T & x, Eigen::Index K, return_type_t<T> & lp);

template <bool Jacobian, typename T, require_std_vector_t<T> *  = nullptr>
 inline  auto cov_matrix_constrain_lkj(const T & x, size_t k, return_type_t<T> & lp);

template <typename T, require_std_vector_t<T> *  = nullptr>
  auto cov_matrix_free(const T & x);

template <typename T, require_std_vector_t<T> *  = nullptr>
  auto cov_matrix_free_lkj(const T & x);

template <bool Jacobian, typename T, require_std_vector_t<T> *  = nullptr>
 inline  auto ordered_constrain(const T & x, return_type_t<T> & lp);

template <typename T, require_std_vector_t<T> *  = nullptr>
  auto ordered_free(const T & x);

template <bool Jacobian, typename T, require_std_vector_t<T> *  = nullptr>
 inline  auto positive_ordered_constrain(const T & x, return_type_t<T> & lp);

template <typename T, require_std_vector_t<T> *  = nullptr>
  auto positive_ordered_free(const T & x);

template <bool Jacobian, typename T, require_std_vector_t<T> *  = nullptr>
 inline  auto simplex_constrain(const T & y, return_type_t<T> & lp);

template <typename T, require_std_vector_t<T> *  = nullptr>
  auto simplex_free(const T & x);

template <bool Jacobian, typename T, require_std_vector_t<T> *  = nullptr>
 inline  auto stochastic_column_constrain(const T & y, return_type_t<T> & lp);

template <typename T, require_std_vector_t<T> *  = nullptr>
 inline  auto stochastic_column_free(const T & y);

template <bool Jacobian, typename T, require_std_vector_t<T> *  = nullptr>
 inline  auto stochastic_row_constrain(const T & y, return_type_t<T> & lp);

template <typename T, require_std_vector_t<T> *  = nullptr>
 inline  auto stochastic_row_free(const T & y);

template <typename T, require_std_vector_t<T> *  = nullptr>
  auto sum_to_zero_free(const T & z);

template <bool Jacobian, typename T, require_std_vector_t<T> *  = nullptr>
 inline  auto unit_vector_constrain(const T & y, return_type_t<T> & lp);

template <typename T, require_std_vector_t<T> *  = nullptr>
  auto unit_vector_free(const T & x);

template <typename T, require_std_vector_st<stan::is_var, T> *  = nullptr>
 inline  auto log_softmax(const T & x);

template <typename T, require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T> *  = nullptr, require_not_var_matrix_t<T> *  = nullptr>
 inline  auto std_normal_log_qf(const T & x);

template <typename T, require_vector_st<stan::is_fvar, T> *  = nullptr>
 inline  auto log_softmax(const T & x);

// Functions Parsed: 2405

} // namespace math
} // namespace stan
#endif // STAN_MATH_VECTORIZE_VECTORIZE_DECL_HPP

