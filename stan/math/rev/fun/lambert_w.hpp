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

/**
 * Internal helper struct for Lambert W function on W0 branch.
 */
class lambertw0_vari : public op_v_vari {
 public:
  explicit lambertw0_vari(vari* avi) : op_v_vari(lambert_w0(avi->val_), avi) {}
  void chain() { avi_->adj_ += (adj_ / (avi_->val_ + exp(val_))); }
};

/**
 * Internal helper struct for Lambert W function on W-1 branch.
 */
class lambertwm1_vari : public op_v_vari {
 public:
  explicit lambertwm1_vari(vari* avi)
      : op_v_vari(lambert_wm1(avi->val_), avi) {}
  void chain() { avi_->adj_ += (adj_ / (avi_->val_ + exp(val_))); }
};
}  // namespace internal

/**
 * Return the Lambert W function on W0 branch applied to the specified variable.
 *
 * @param a Variable argument.
 * @return the Lambert W function (W0 branch) applied to the specified argument.
 */
inline var lambert_w0(const var& a) {
  return var(new internal::lambertw0_vari(a.vi_));
}

/**
 * Return the Lambert W function on W-1 branch applied to the specified
 * variable.
 *
 * @param a Variable argument.
 * @return the Lambert W function (W-1 branch) applied to the specified
 * argument.
 */
inline var lambert_wm1(const var& a) {
  return var(new internal::lambertwm1_vari(a.vi_));
}

}  // namespace math
}  // namespace stan

#endif
