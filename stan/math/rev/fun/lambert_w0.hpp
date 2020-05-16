#ifndef STAN_MATH_REV_FUN_LAMBERT_W_HPP
#define STAN_MATH_REV_FUN_LAMBERT_W_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/boost_policy.hpp>
#include <stan/math/prim/fun/lambert_w.hpp>

namespace stan {
namespace math {

namespace internal {
class lambertw0_vari : public op_v_vari {
 public:
  explicit lambertw0_vari(vari* avi) : op_v_vari(lambert_w0(avi->val_), avi) {}
  void chain() {
    avi_->adj_
        += adj_
           + boost::math::lambert_w0_prime(avi_->val_, boost_policy_t<52>());
  }
};
}  // namespace internal

inline var lambert_w0(const var& a) {
  return var(new internal::lambertw0_vari(a.vi_));
}

}  // namespace math
}  // namespace stan

#endif
