#ifndef STAN_MATH_REV_FUN_INC_BETA_HPP
#define STAN_MATH_REV_FUN_INC_BETA_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/beta.hpp>
#include <stan/math/rev/fun/digamma.hpp>
#include <stan/math/rev/fun/pow.hpp>
#include <stan/math/prim/fun/grad_reg_inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {

class inc_beta_vvv_vari : public op_vvv_vari {
 public:
  inc_beta_vvv_vari(vari* avi, vari* bvi, vari* cvi)
      : op_vvv_vari(inc_beta(avi->val_, bvi->val_, cvi->val_), avi, bvi, cvi) {}
  void chain() {
    double d_a;
    double d_b;
    const double beta_ab = beta(avi_->val_, bvi_->val_);
    grad_reg_inc_beta(d_a, d_b, avi_->val_, bvi_->val_, cvi_->val_,
                      digamma(avi_->val_), digamma(bvi_->val_),
                      digamma(avi_->val_ + bvi_->val_), beta_ab);

    avi_->adj_ += adj_ * d_a;
    bvi_->adj_ += adj_ * d_b;
    cvi_->adj_ += adj_ * std::pow(1 - cvi_->val_, bvi_->val_ - 1)
                  * std::pow(cvi_->val_, avi_->val_ - 1) / beta_ab;
  }
};

}  // namespace internal

inline var inc_beta(const var& a, const var& b, const var& c) {
  return var(new internal::inc_beta_vvv_vari(a.vi_, b.vi_, c.vi_));
}

/*
template <typename Ta, typename Tb, typename Tc,
          require_all_stan_scalar_t<Ta, Tb, Tc>* = nullptr,
          require_any_var_t<Ta, Tb, Tc>* = nullptr>
inline var inc_beta(const Ta& a, const Tb& b, const Tc& c) {
  double a_dbl = value_of(a);
  double b_dbl = value_of(b);
  double c_dbl = value_of(c);
  double rtn = inc_beta(a_dbl, b_dbl, c_dbl);
  auto grad_tuple = grad_2F1(a+b, 1, a+1, c_dbl);
  return make_callback_var(rtn, [a, b, c, a_dbl, b_dbl, c_dbl, rtn, grad_tuple](auto& vi) mutable {
    if (!is_constant<Ta>::value) {
      double da =
        rtn / a_dbl
        + rtn * log(c_dbl)
        - rtn * (digamma(a_dbl) - digamma(a_dbl + b_dbl))
        + rtn * (std::get<0>(grad_tuple) + std::get<2>(grad_tuple));
      forward_as<var>(a).adj() += vi.adj() * da;
    }
    if (!is_constant<Tb>::value) {
      double db =
        rtn * log1m(c_dbl)
        - rtn * (digamma(b_dbl) - digamma(a_dbl + b_dbl))
        + rtn * std::get<0>(grad_tuple);
      forward_as<var>(b).adj() += vi.adj() * db;
    }
    if (!is_constant<Tc>::value) {
      double beta_ab = beta(a_dbl, b_dbl);
      forward_as<var>(c).adj() += vi.adj() * rtn * std::pow(1 - c_dbl, b_dbl - 1) * std::pow(c_dbl, a_dbl - 1) / beta_ab;
      //forward_as<var>(c).adj() += vi.adj() * std::get<2>(grad_tuple);
    }
  });
}
*/

}  // namespace math
}  // namespace stan
#endif
