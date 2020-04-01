#ifndef STAN_MATH_REV_FUN_CONV_GAUS_LINE
#define STAN_MATH_REV_FUN_CONV_GAUS_LINE

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/conv_gaus_line.hpp>

namespace stan {
namespace math {

namespace internal {
class conv_gaus_line_vari : public op_v_vari {
  double t0_;
  double t1_;
  double a_;
  double b_;
  double sig2_;

 public:
  explicit conv_gaus_line_vari(double t0, double t1, double a, double b,
                               vari* avi, double sig2)
      : t0_(t0),
        t1_(t1),
        a_(a),
        b_(b),
        sig2_(sig2),
        op_v_vari(conv_gaus_line(t0, t1, a, b, avi->val_, sig2), avi) {}

  void chain() {
    avi_->adj_
        += adj_ * der_conv_gaus_line(t0_, t1_, a_, b_, avi_->val_, sig2_);
  }
};

}  // namespace internal

inline var conv_gaus_line(double t0, double t1, double a, double b,
                          const var& x, double sig2) {
  return var(new internal::conv_gaus_line_vari(t0, t1, a, b, x.vi_, sig2));
}

}  // namespace math
}  // namespace stan

#endif
