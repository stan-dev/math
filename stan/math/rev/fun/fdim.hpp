#ifndef STAN_MATH_REV_FUN_FDIM_HPP
#define STAN_MATH_REV_FUN_FDIM_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>

namespace stan {
namespace math {

namespace internal {
class fdim_vv_vari : public op_vv_vari {
 public:
  fdim_vv_vari(vari* avi, vari* bvi)
      : op_vv_vari(avi->val_ - bvi->val_, avi, bvi) {}
  void chain() {
    if (unlikely(is_any_nan(avi_->val_, bvi_->val_))) {
      avi_->adj_ = NOT_A_NUMBER;
      bvi_->adj_ = NOT_A_NUMBER;
    } else {
      avi_->adj_ += adj_;
      bvi_->adj_ -= adj_;
    }
  }
};

class fdim_vd_vari : public op_vd_vari {
 public:
  fdim_vd_vari(vari* avi, double b) : op_vd_vari(avi->val_ - b, avi, b) {}
  void chain() {
    if (unlikely(is_any_nan(avi_->val_, bd_))) {
      avi_->adj_ = NOT_A_NUMBER;
    } else {
      avi_->adj_ += adj_;
    }
  }
};

class fdim_dv_vari : public op_dv_vari {
 public:
  fdim_dv_vari(double a, vari* bvi) : op_dv_vari(a - bvi->val_, a, bvi) {}
  void chain() {
    if (unlikely(is_any_nan(bvi_->val_, ad_))) {
      bvi_->adj_ = NOT_A_NUMBER;
    } else {
      bvi_->adj_ -= adj_;
    }
  }
};
}  // namespace internal

/**
 * Return the positive difference between the first variable's the value
 * and the second's (C99, C++11).
 *
 * The function values and derivatives are defined by
 *
   \f[
   \mbox{fdim}(x, y) =
   \begin{cases}
     x-y & \mbox{if } x  > y \\[6pt]
     0 & \mbox{otherwise} \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{fdim}(x, y)}{\partial x} =
   \begin{cases}
     1 & \mbox{if } x > y \\[6pt]
     0 & \mbox{otherwise} \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{fdim}(x, y)}{\partial y} =
   \begin{cases}
    -1 & \mbox{if } x > y \\[6pt]
     0 & \mbox{otherwise} \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a First variable.
 * @param b Second variable.
 * @return The positive difference between the first and second
 * variable.
 */
inline var fdim(const var& a, const var& b) {
  // reversed test to get NaN vals automatically in second case
  return (a.vi_->val_ <= b.vi_->val_)
             ? var(new vari(0.0))
             : var(new internal::fdim_vv_vari(a.vi_, b.vi_));
}

/**
 * Return the positive difference between the first value and the
 * value of the second variable (C99, C++11).
 *
 * See <code>fdim(var, var)</code> for definitions of values and
 * derivatives.
 *
 * @param a First value.
 * @param b Second variable.
 * @return The positive difference between the first and second
 * arguments.
 */
inline var fdim(double a, const var& b) {
  // reversed test to get NaN vals automatically in second case
  return a <= b.vi_->val_ ? var(new vari(0.0))
                          : var(new internal::fdim_dv_vari(a, b.vi_));
}

/**
 * Return the positive difference between the first variable's value
 * and the second value (C99, C++11).
 *
 * See <code>fdim(var, var)</code> for definitions of values and
 * derivatives.
 *
 * @param a First value.
 * @param b Second variable.
 * @return The positive difference between the first and second arguments.
 */
inline var fdim(const var& a, double b) {
  // reversed test to get NaN vals automatically in second case
  return a.vi_->val_ <= b ? var(new vari(0.0))
                          : var(new internal::fdim_vd_vari(a.vi_, b));
}

}  // namespace math
}  // namespace stan
#endif
