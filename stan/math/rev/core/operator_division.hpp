#ifndef STAN_MATH_REV_CORE_OPERATOR_DIVISION_HPP
#define STAN_MATH_REV_CORE_OPERATOR_DIVISION_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/core/operator_division.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/meta/is_vari.hpp>
#include <stan/math/rev/core/std_complex.hpp>
#include <stan/math/rev/core/op_vari.hpp>
#include <stan/math/rev/core/operator_addition.hpp>
#include <stan/math/rev/core/operator_multiplication.hpp>
#include <stan/math/rev/core/operator_subtraction.hpp>
#include <stan/math/rev/meta/apply_scalar_unary.hpp>
#include <complex>
#include <type_traits>

namespace stan {
namespace math {

namespace internal {
template <typename VariVal, typename Vari1, typename Vari2, typename = void>
class divide_vari {};

// (dividend/divisor)' = dividend' * (1 / divisor) - divisor' * (dividend /
// [divisor * divisor])
template <typename VariVal, typename Vari1, typename Vari2>
class divide_vari<VariVal, Vari1, Vari2, require_all_vari_t<Vari1, Vari2>>
    : public op_vari<VariVal, Vari1*, Vari2*> {
  using op_vari<VariVal, Vari1*, Vari2*>::avi;
  using op_vari<VariVal, Vari1*, Vari2*>::bvi;

 public:
  divide_vari(Vari1* dividend_vi, Vari2* divisor_vi)
      : op_vari<VariVal, Vari1*, Vari2*>(dividend_vi->val_ / divisor_vi->val_,
                                         dividend_vi, divisor_vi) {}
  void chain() {
    if (unlikely(is_any_nan(avi()->val_, bvi()->val_))) {
      avi()->adj_ = NOT_A_NUMBER;
      bvi()->adj_ = NOT_A_NUMBER;
    } else {
      avi()->adj_ += this->adj_ / bvi()->val_;
      bvi()->adj_ -= this->adj_ * avi()->val_ / (bvi()->val_ * bvi()->val_);
    }
  }
};

template <typename VariVal, typename Vari, typename Arith>
class divide_vari<
    VariVal, Vari, Arith,
    require_t<conjunction<std::is_arithmetic<Arith>, is_vari<Vari>>>>
    : public op_vari<VariVal, Vari*, Arith> {
  using op_vari<VariVal, Vari*, Arith>::avi;
  using op_vari<VariVal, Vari*, Arith>::bd;

 public:
  divide_vari(Vari* dividend_vi, Arith divisor)
      : op_vari<VariVal, Vari*, Arith>(dividend_vi->val_ / divisor, dividend_vi,
                                       divisor) {}
  void chain() {
    if (unlikely(is_any_nan(avi()->val_, bd()))) {
      avi()->adj_ = NOT_A_NUMBER;
    } else {
      avi()->adj_ += this->adj_ / bd();
    }
  }
};

template <typename VariVal, typename Arith, typename Vari>
class divide_vari<
    VariVal, Arith, Vari,
    require_t<conjunction<std::is_arithmetic<Arith>, is_vari<Vari>>>>
    : public op_vari<VariVal, Arith, Vari*> {
  using op_vari<VariVal, Arith, Vari*>::ad;
  using op_vari<VariVal, Arith, Vari*>::bvi;

 public:
  divide_vari(const Arith& dividend, Vari* divisor_vi)
      : op_vari<VariVal, Arith, Vari*>(dividend / divisor_vi->val_, dividend,
                                       divisor_vi) {}
  void chain() {
    bvi()->adj_ -= this->adj_ * ad() / (bvi()->val_ * bvi()->val_);
  }
};
}  // namespace internal

/**
 * Division operator for two variables (C++).
 *
 * The partial derivatives for the variables are
 *
 * \f$\frac{\partial}{\partial x} (x/y) = 1/y\f$, and
 *
 * \f$\frac{\partial}{\partial y} (x/y) = -x / y^2\f$.
 *
   \f[
   \mbox{operator/}(x, y) =
   \begin{cases}
     \frac{x}{y} & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
   \f[
   \frac{\partial\, \mbox{operator/}(x, y)}{\partial x} =
   \begin{cases}
     \frac{1}{y} & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
   \f[
   \frac{\partial\, \mbox{operator/}(x, y)}{\partial y} =
   \begin{cases}
     -\frac{x}{y^2} & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param dividend First variable operand.
 * @param divisor Second variable operand.
 * @return Variable result of dividing the first variable by the
 * second.
 */
template <typename T>
inline var_value<T> operator/(const var_value<T>& dividend,
                              const var_value<T>& divisor) {
  return {new internal::divide_vari<T, vari_value<T>, vari_value<T>>(
      dividend.vi_, divisor.vi_)};
}

/**
 * Division operator for dividing a variable by a scalar (C++).
 *
 * The derivative with respect to the variable is
 *
 * \f$\frac{\partial}{\partial x} (x/c) = 1/c\f$.
 *
 * @tparam Var value type of a var
 * @tparam Arith An arithmetic type
 * @param dividend Variable operand.
 * @param divisor Scalar operand.
 * @return Variable result of dividing the variable by the scalar.
 */
template <typename T, typename Arith, require_arithmetic_t<Arith>...>
inline var_value<T> operator/(const var_value<T>& dividend,
                              const Arith& divisor) {
  if (divisor == 1.0) {
    return dividend;
  }
  return {new internal::divide_vari<T, vari_value<T>, Arith>(dividend.vi_,
                                                             divisor)};
}

/**
 * Division operator for dividing a scalar by a variable (C++).
 *
 * The derivative with respect to the variable is
 *
 * \f$\frac{d}{d y} (c/y) = -c / y^2\f$.
 *
 * @tparam Arith An arithmetic type
 * @param dividend Scalar operand.
 * @param divisor Variable operand.
 * @return Quotient of the dividend and divisor.
 */
template <typename T, typename Arith, require_arithmetic_t<Arith>...>
inline var_value<T> operator/(const Arith& dividend,
                              const var_value<T>& divisor) {
  return {new internal::divide_vari<T, Arith, vari_value<T>>(dividend,
                                                             divisor.vi_)};
}

inline std::complex<var> operator/(const std::complex<var>& x1,
                                   const std::complex<var>& x2) {
  return internal::complex_divide(x1, x2);
}

}  // namespace math
}  // namespace stan
#endif
