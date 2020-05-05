#ifndef STAN_MATH_REV_CORE_OPERATOR_MULTIPLICATION_HPP
#define STAN_MATH_REV_CORE_OPERATOR_MULTIPLICATION_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/vv_vari.hpp>
#include <stan/math/rev/core/vd_vari.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <stan/math/prim/fun/isinf.hpp>
#include <stan/math/prim/fun/isnan.hpp>

namespace stan {
namespace math {

namespace internal {
template <typename VariVal, typename Vari1, typename Vari2, typename = void>
class multiply_vari;

template <typename VariVal, typename Vari1, typename Vari2>
class multiply_vari<VariVal, Vari1, Vari2, require_all_vari_t<Vari1, Vari2>> :
  public op_vari<VariVal, Vari1*, Vari2*> {
 public:
  multiply_vari(Vari1* avi, Vari2* bvi)
      : op_vari<VariVal, Vari1*, Vari2*>(avi->val_ * bvi->val_, avi, bvi) {}
  void chain() {
  if (likely(is_not_nan(std::get<0>(this->vi())->val_) &&
                      is_not_nan(std::get<1>(this->vi())->val_))) {
      this->avi()->adj_ += this->bvi()->val_ * this->adj_;
      this->bvi()->adj_ += this->avi()->val_ * this->adj_;
    } else {
      fill(std::get<0>(this->vi())->adj_, NOT_A_NUMBER);
      fill(std::get<1>(this->vi())->adj_, NOT_A_NUMBER);
    }
  }
};

template <typename VariVal, typename Vari, typename Arith>
class multiply_vari<VariVal, Vari, Arith, require_vt_arithmetic<Arith>> : public op_vari<VariVal, Vari*, Arith> {
 public:
  multiply_vari(Vari* avi, Arith b) : op_vari<VariVal, Vari*, Arith>(avi->val_ * b, avi, b) {}
  void chain() {
    if (unlikely(is_any_nan(this->avi()->val_, this->bvi()))) {
      fill(this->avi()->adj_, NOT_A_NUMBER);
    } else {
      this->avi()->adj_ += this->adj_ * this->bvi();
    }
  }
};
}  // namespace internal

/**
 * Multiplication operator for two variables (C++).
 *
 * The partial derivatives are
 *
 * \f$\frac{\partial}{\partial x} (x * y) = y\f$, and
 *
 * \f$\frac{\partial}{\partial y} (x * y) = x\f$.
 *
   \f[
   \mbox{operator*}(x, y) =
   \begin{cases}
     xy & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{operator*}(x, y)}{\partial x} =
   \begin{cases}
     y & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{operator*}(x, y)}{\partial y} =
   \begin{cases}
     x & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of multiplying operands.
 */
template <typename T>
inline var_type<T> operator*(const var_type<T>& a, const var_type<T>& b) {
  return {new internal::multiply_vari<T, vari_type<T>, vari_type<T>>(a.vi_, b.vi_)};
}

// Just shoving this here for now
namespace internal {
  template <typename T1, typename T2, require_all_arithmetic_t<T1, T2>* = nullptr>
  bool is_any_equal(T1 x, T2 y) {
    return x == y;
  }
  template <typename T1, typename T2, require_eigen_t<T1>* = nullptr, require_arithmetic_t<T2>* = nullptr>
  bool is_any_equal(const T1& x, T2 y) {
    return (x.array() == y).any();
  }

}
/**
 * Multiplication operator for a variable and a scalar (C++).
 *
 * The partial derivative for the variable is
 *
 * \f$\frac{\partial}{\partial x} (x * c) = c\f$, and
 *
 * @tparam Arith An arithmetic type
 * @param a Variable operand.
 * @param b Scalar operand.
 * @return Variable result of multiplying operands.
 */
template <typename T, typename Arith, require_vt_arithmetic<Arith>...>
inline var_type<T> operator*(const var_type<T>& a, Arith b) {
  if (internal::is_any_equal(b, 1.0)) {
    return a;
  }
  return {new internal::multiply_vari<T, vari_type<T>, Arith>(a.vi_, b)};
}

/**
 * Multiplication operator for a scalar and a variable (C++).
 *
 * The partial derivative for the variable is
 *
 * \f$\frac{\partial}{\partial y} (c * y) = c\f$.
 *
 * @tparam Arith An arithmetic type
 * @param a Scalar operand.
 * @param b Variable operand.
 * @return Variable result of multiplying the operands.
 */


template <typename T, typename Arith, require_vt_arithmetic<Arith>...>
inline var_type<T> operator*(Arith a, const var_type<T>& b) {
  if (internal::is_any_equal(a, 1.0)) {
    return b;
  }
  return {new internal::multiply_vari<double, vari, Arith>(b.vi_, a)};  // by symmetry
}

}  // namespace math
}  // namespace stan
#endif
