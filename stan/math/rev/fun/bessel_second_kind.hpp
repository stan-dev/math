#ifndef STAN_MATH_REV_FUN_BESSEL_SECOND_KIND_HPP
#define STAN_MATH_REV_FUN_BESSEL_SECOND_KIND_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/bessel_second_kind.hpp>

namespace stan {
namespace math {

namespace internal {

class bessel_second_kind_dv_vari : public op_dv_vari {
 public:
  bessel_second_kind_dv_vari(int a, vari* bvi)
      : op_dv_vari(bessel_second_kind(a, bvi->val_), a, bvi) {}
  void chain() {
    bvi_->adj_ += adj_
                  * (ad_ * bessel_second_kind(ad_, bvi_->val_) / bvi_->val_
                     - bessel_second_kind(ad_ + 1, bvi_->val_));
  }
};
}  // namespace internal

inline var bessel_second_kind(int v, const var& a) {
  return make_callback_var(bessel_second_kind(v, a.val()),
    [v, a](auto& vi) mutable {
      a.adj() += vi.adj()
                    * (v * bessel_second_kind(v, a.val()) / a.val()
                       - bessel_second_kind(v + 1, a.val()));
    });
}

}  // namespace math
}  // namespace stan
#endif
