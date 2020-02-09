#ifndef STAN_MATH_REV_CORE_OPERATOR_ADDITION_HPP
#define STAN_MATH_REV_CORE_OPERATOR_ADDITION_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/op_vari.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <type_traits>

namespace stan {
namespace math {

namespace internal {

template <typename T>
struct is_vari : bool_constant<std::is_base_of<
                     vari, std::remove_pointer_t<std::decay_t<T>>>::value> {};

template <typename T>
using require_vari_t = require_t<is_vari<T>>;

template <typename... Types>
using require_all_vari_t = require_all_t<is_vari<Types>...>;

template <typename T1, typename T2, typename = void>
class add_vari {
  static_assert(1, "If you see this please report a bug!");
};

template <typename T1, typename T2>
class add_vari<T1, T2, require_all_vari_t<T1, T2>> : public op_vari<T1, T2> {
 public:
  add_vari(T1 avi, T2 bvi) : op_vari<T2, T1>(avi->val_ + bvi->val_, avi, bvi) {}
  void chain() {
    if (unlikely(is_any_nan(std::get<0>(this->vi())->val_,
                            std::get<1>(this->vi())->val_))) {
      std::get<0>(this->vi())->adj_ = NOT_A_NUMBER;
      std::get<1>(this->vi())->adj_ = NOT_A_NUMBER;
    } else {
      std::get<0>(this->vi())->adj_ += this->adj_;
      std::get<1>(this->vi())->adj_ += this->adj_;
    }
  }
};

template <typename T1, typename T2>
class add_vari<T1, T2,
               require_t<conjunction<is_vari<T1>, std::is_floating_point<T2>>>>
    : public op_vari<T1, T2> {
 public:
  add_vari(T1 avi, T2 b) : op_vari<T1, T2>(avi->val_ + b, avi, b) {}
  void chain() {
    if (unlikely(is_any_nan(std::get<0>(this->vi())->val_,
                            std::get<1>(this->vi())))) {
      std::get<0>(this->vi())->adj_ = NOT_A_NUMBER;
    } else {
      std::get<0>(this->vi())->adj_ += this->adj_;
    }
  }
};
}  // namespace internal

/**
 * Addition operator for variables (C++).
 *
 * The partial derivatives are defined by
 *
 * \f$\frac{\partial}{\partial x} (x+y) = 1\f$, and
 *
 * \f$\frac{\partial}{\partial y} (x+y) = 1\f$.
 *
 *
   \f[
   \mbox{operator+}(x, y) =
   \begin{cases}
     x+y & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{operator+}(x, y)}{\partial x} =
   \begin{cases}
     1 & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{operator+}(x, y)}{\partial y} =
   \begin{cases}
     1 & \mbox{if } -\infty\leq x, y \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a First variable operand.
 * @param b Second variable operand.
 * @return Variable result of adding two variables.
 */
inline var operator+(var a, var b) {
  return {
      new internal::add_vari<decltype(a.vi_), decltype(b.vi_)>(a.vi_, b.vi_)};
}

/**
 * Addition operator for variable and scalar (C++).
 *
 * The derivative with respect to the variable is
 *
 * \f$\frac{d}{dx} (x + c) = 1\f$.
 *
 * @tparam Arith An arithmetic type
 * @param a First variable operand.
 * @param b Second scalar operand.
 * @return Result of adding variable and scalar.
 */
template <typename Arith, require_arithmetic_t<Arith>...>
inline var operator+(var a, Arith b) {
  if (b == 0.0) {
    return a;
  }
  return {new internal::add_vari<decltype(a.vi_), double>(a.vi_, b)};
}

/**
 * Addition operator for scalar and variable (C++).
 *
 * The derivative with respect to the variable is
 *
 * \f$\frac{d}{dy} (c + y) = 1\f$.
 *
 * @tparam Arith An arithmetic type
 * @param a First scalar operand.
 * @param b Second variable operand.
 * @return Result of adding variable and scalar.
 */
template <typename Arith, require_arithmetic_t<Arith>...>
inline var operator+(Arith a, var b) {
  if (a == 0.0) {
    return b;
  }
  return {new internal::add_vari<decltype(b.vi_), double>(b.vi_,
                                                          a)};  // by symmetry
}

}  // namespace math
}  // namespace stan
#endif
