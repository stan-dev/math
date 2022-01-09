#ifndef STAN_MATH_FWD_PROB_EXPONENTIAL_QF_HPP
#define STAN_MATH_FWD_PROB_EXPONENTIAL_QF_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/prob/exponential_qf.hpp>

namespace stan {
namespace math {

template <typename Tp, typename Tbeta,
          require_any_st_fvar<Tp, Tbeta>* = nullptr>
inline auto exponential_qf(const Tp& p, const Tbeta& beta) {
  using vector_t = decltype(exponential_qf(value_of(p),
                                           value_of(beta)));
  using fvar_t = return_type_t<Tp, Tbeta>;

  auto p_val = value_of(p);
  auto beta_val = value_of(beta);
  promote_scalar_t<fvar_t, vector_t>
    ret(exponential_qf(p_val, beta_val));
  ret.d();// = vector_t(0.0);
  if(is_fvar<Tp>::value) {
    as_array_or_scalar(ret.d()) +=
      as_array_or_scalar(forward_as<promote_scalar_t<fvar_t, Tp>>(p).d())
      * inv(as_array_or_scalar(beta_val)
            - as_array_or_scalar(beta_val)
            * as_array_or_scalar(p_val));
  }
  if(is_fvar<Tbeta>::value) {
    as_array_or_scalar(ret.d()) +=
      as_array_or_scalar(forward_as<promote_scalar_t<fvar_t, Tbeta>>(beta).d())
      * log1m(as_array_or_scalar(p_val)) * as_array_or_scalar(beta_val);
  }

  return ret;
}


}  // namespace math
}  // namespace stan
#endif
