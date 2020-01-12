#ifndef STAN_MATH_REV_FUN_FMA_HPP
#define STAN_MATH_REV_FUN_FMA_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/fma.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>

namespace stan {
namespace math {

namespace internal {
class fma_vvv_vari : public op_vvv_vari {
 public:
  fma_vvv_vari(vari* avi, vari* bvi, vari* cvi)
      : op_vvv_vari(fma(avi->val_, bvi->val_, cvi->val_), avi, bvi, cvi) {}
  void chain() {
    if (unlikely(is_any_nan(avi_->val_, bvi_->val_, cvi_->val_))) {
      avi_->adj_ = NOT_A_NUMBER;
      bvi_->adj_ = NOT_A_NUMBER;
      cvi_->adj_ = NOT_A_NUMBER;
    } else {
      avi_->adj_ += adj_ * bvi_->val_;
      bvi_->adj_ += adj_ * avi_->val_;
      cvi_->adj_ += adj_;
    }
  }
};

class fma_vvd_vari : public op_vvd_vari {
 public:
  fma_vvd_vari(vari* avi, vari* bvi, double c)
      : op_vvd_vari(fma(avi->val_, bvi->val_, c), avi, bvi, c) {}
  void chain() {
    if (unlikely(is_any_nan(avi_->val_, bvi_->val_, cd_))) {
      avi_->adj_ = NOT_A_NUMBER;
      bvi_->adj_ = NOT_A_NUMBER;
    } else {
      avi_->adj_ += adj_ * bvi_->val_;
      bvi_->adj_ += adj_ * avi_->val_;
    }
  }
};

class fma_vdv_vari : public op_vdv_vari {
 public:
  fma_vdv_vari(vari* avi, double b, vari* cvi)
      : op_vdv_vari(fma(avi->val_, b, cvi->val_), avi, b, cvi) {}
  void chain() {
    if (unlikely(is_any_nan(avi_->val_, cvi_->val_, bd_))) {
      avi_->adj_ = NOT_A_NUMBER;
      cvi_->adj_ = NOT_A_NUMBER;
    } else {
      avi_->adj_ += adj_ * bd_;
      cvi_->adj_ += adj_;
    }
  }
};

class fma_vdd_vari : public op_vdd_vari {
 public:
  fma_vdd_vari(vari* avi, double b, double c)
      : op_vdd_vari(fma(avi->val_, b, c), avi, b, c) {}
  void chain() {
    if (unlikely(is_any_nan(avi_->val_, bd_, cd_))) {
      avi_->adj_ = NOT_A_NUMBER;
    } else {
      avi_->adj_ += adj_ * bd_;
    }
  }
};

class fma_ddv_vari : public op_ddv_vari {
 public:
  fma_ddv_vari(double a, double b, vari* cvi)
      : op_ddv_vari(fma(a, b, cvi->val_), a, b, cvi) {}
  void chain() {
    if (unlikely(is_any_nan(cvi_->val_, ad_, bd_))) {
      cvi_->adj_ = NOT_A_NUMBER;
    } else {
      cvi_->adj_ += adj_;
    }
  }
};
}  // namespace internal

/**
 * The fused multiply-add function for three variables (C99).
 * This function returns the product of the first two arguments
 * plus the third argument.
 *
 * The partial derivatives are
 *
 * \f$\frac{\partial}{\partial x} (x * y) + z = y\f$, and
 *
 * \f$\frac{\partial}{\partial y} (x * y) + z = x\f$, and
 *
 * \f$\frac{\partial}{\partial z} (x * y) + z = 1\f$.
 *
 * @param a First multiplicand.
 * @param b Second multiplicand.
 * @param c Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
inline var fma(const var& a, const var& b, const var& c) {
  return var(new internal::fma_vvv_vari(a.vi_, b.vi_, c.vi_));
}

/**
 * The fused multiply-add function for two variables and a value
 * (C99).  This function returns the product of the first two
 * arguments plus the third argument.
 *
 * The partial derivatives are
 *
 * \f$\frac{\partial}{\partial x} (x * y) + c = y\f$, and
 *
 * \f$\frac{\partial}{\partial y} (x * y) + c = x\f$.
 *
 * @tparam Tc type of the summand
 * @param a First multiplicand.
 * @param b Second multiplicand.
 * @param c Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Tc, typename = require_arithmetic_t<Tc>>
inline var fma(const var& a, const var& b, Tc&& c) {
  return var(new internal::fma_vvd_vari(a.vi_, b.vi_, c));
}

/**
 * The fused multiply-add function for a variable, value, and
 * variable (C99).  This function returns the product of the first
 * two arguments plus the third argument.
 *
 * The partial derivatives are
 *
 * \f$\frac{\partial}{\partial x} (x * c) + z = c\f$, and
 *
 * \f$\frac{\partial}{\partial z} (x * c) + z = 1\f$.
 *
 * @tparam Ta type of the first multiplicand
 * @tparam Tb type of the second multiplicand
 * @tparam Tc type of the summand
 *
 * @param a First multiplicand.
 * @param b Second multiplicand.
 * @param c Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Ta, typename Tb, typename Tc,
          typename = require_arithmetic_t<Tb>,
          typename = require_all_var_t<Ta, Tc>>
inline var fma(Ta&& a, Tb&& b, Tc&& c) {
  return var(new internal::fma_vdv_vari(a.vi_, b, c.vi_));
}

/**
 * The fused multiply-add function for a variable and two values
 * (C99).  This function returns the product of the first two
 * arguments plus the third argument.
 *
 * The double-based version
 * <code>::%fma(double, double, double)</code> is defined in
 * <code>&lt;cmath&gt;</code>.
 *
 * The derivative is
 *
 * \f$\frac{d}{d x} (x * c) + d = c\f$.
 *
 * @tparam Tb type of the second multiplicand
 * @tparam Tc type of the summand
 *
 * @param a First multiplicand.
 * @param b Second multiplicand.
 * @param c Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Tb, typename Tc, typename = require_all_arithmetic_t<Tb, Tc>>
inline var fma(const var& a, Tb&& b, Tc&& c) {
  return var(new internal::fma_vdd_vari(a.vi_, b, c));
}

/**
 * The fused multiply-add function for a value, variable, and
 * value (C99).  This function returns the product of the first
 * two arguments plus the third argument.
 *
 * The derivative is
 *
 * \f$\frac{d}{d y} (c * y) + d = c\f$, and
 *
 * @tparam Ta type of the first multiplicand
 * @tparam Tc type of the summand
 *
 * @param a First multiplicand.
 * @param b Second multiplicand.
 * @param c Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Ta, typename Tc, typename = require_all_arithmetic_t<Ta, Tc>>
inline var fma(Ta&& a, const var& b, Tc&& c) {
  return var(new internal::fma_vdd_vari(b.vi_, a, c));
}

/**
 * The fused multiply-add function for two values and a variable,
 * and value (C99).  This function returns the product of the
 * first two arguments plus the third argument.
 *
 * The derivative is
 *
 * \f$\frac{\partial}{\partial z} (c * d) + z = 1\f$.
 *
 * @tparam Ta type of the first multiplicand
 * @tparam Tb type of the second multiplicand
 *
 * @param a First multiplicand.
 * @param b Second multiplicand.
 * @param c Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Ta, typename Tb, typename = require_all_arithmetic_t<Ta, Tb>>
inline var fma(Ta&& a, Tb&& b, const var& c) {
  return var(new internal::fma_ddv_vari(a, b, c.vi_));
}

/**
 * The fused multiply-add function for a value and two variables
 * (C99).  This function returns the product of the first two
 * arguments plus the third argument.
 *
 * The partial derivaties are
 *
 * \f$\frac{\partial}{\partial y} (c * y) + z = c\f$, and
 *
 * \f$\frac{\partial}{\partial z} (c * y) + z = 1\f$.
 *
 * @tparam Ta type of the first multiplicand
 * @param a First multiplicand.
 * @param b Second multiplicand.
 * @param c Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename Ta, typename = require_arithmetic_t<Ta>>
inline var fma(Ta&& a, const var& b, const var& c) {
  return var(new internal::fma_vdv_vari(b.vi_, a, c.vi_));  // a-b symmetry
}

}  // namespace math
}  // namespace stan
#endif
