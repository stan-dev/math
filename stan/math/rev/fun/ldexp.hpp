#ifndef STAN_MATH_REV_FUN_LDEXP_HPP
#define STAN_MATH_REV_FUN_LDEXP_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/ldexp.hpp>

namespace stan {
namespace math {

namespace {
class ldexp_vari : public op_vd_vari {
 public:
  explicit ldexp_vari(vari* avi, int b)
      : op_vd_vari(ldexp(avi->val_, b), avi, b) {}
  void chain() { avi_->adj_ += ldexp(adj_, bd_); }
};
}  // namespace

inline var ldexp(const var& a, int b) { return var(new ldexp_vari(a.vi_, b)); }

}  // namespace math
}  // namespace stan
#endif
