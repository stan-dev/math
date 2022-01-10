#ifndef STAN_MATH_FWD_PROB_EXPONENTIAL_QF_HPP
#define STAN_MATH_FWD_PROB_EXPONENTIAL_QF_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/prob/exponential_qf.hpp>

namespace stan {
namespace math {

template <typename Tp, typename Tbeta,
          require_any_st_fvar<Tp, Tbeta>* = nullptr,
          require_all_stan_scalar_t<Tp, Tbeta>* = nullptr>
inline auto exponential_qf(const Tp& p, const Tbeta& beta) {
  using fvar_t = return_type_t<Tp, Tbeta>;

  auto p_val = value_of(p);
  auto beta_val = value_of(beta);
  fvar_t ret(exponential_qf(p_val, beta_val));

  if(!is_constant<Tp>::value) {
    ret.d_ += forward_as<fvar_t>(p).d_
      * inv(beta_val - beta_val * p_val);
  }
  if(!is_constant<Tbeta>::value) {
    ret.d_ += forward_as<fvar_t>(beta).d_
      * log1m(p_val) / square(beta_val);
  }

  return ret;
}

template <typename Tp, typename Tbeta,
          require_any_st_fvar<Tp, Tbeta>* = nullptr,
          require_any_eigen_vector_t<Tp, Tbeta>* = nullptr>
inline auto exponential_qf(const Tp& p, const Tbeta& beta) {
  using vector_t
    = plain_type_t<decltype(exponential_qf(value_of(p), value_of(beta)))>;
  using fvar_t = return_type_t<Tp, Tbeta>;

  auto p_val = value_of(p);
  auto beta_val = value_of(beta);
  vector_t ret = exponential_qf(p_val, beta_val);

  auto p_array = as_array_or_scalar(p_val);
  auto beta_array = as_array_or_scalar(beta_val);
  vector_t d_ = vector_t::Zero(ret.rows(), ret.cols());

  if(!is_constant<Tp>::value) {
    as_array_or_scalar(d_)
      += as_array_or_scalar(forward_as<promote_scalar_t<fvar_t, Tp>>(p).d())
      * inv(beta_array - beta_array * p_array);
  }
  if(!is_constant<Tbeta>::value) {
    as_array_or_scalar(d_)
      += as_array_or_scalar(forward_as<promote_scalar_t<fvar_t,
                                                        Tbeta>>(beta).d())
      * log1m(p_array) / square(beta_array);
  }

  return ret.binaryExpr(d_, [&](const auto& val, const auto& deriv) {
    return fvar_t(val, deriv);
  }).eval();
}

}  // namespace math
}  // namespace stan
#endif
