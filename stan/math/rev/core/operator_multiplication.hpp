#ifndef STAN_MATH_REV_CORE_OPERATOR_MULTIPLICATION_HPP
#define STAN_MATH_REV_CORE_OPERATOR_MULTIPLICATION_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/op_vari.hpp>
#include <stan/math/rev/meta/is_vari.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <stan/math/prim/fun/isinf.hpp>
#include <stan/math/prim/fun/isnan.hpp>

namespace stan {
namespace math {

namespace internal {
/**
 * Base for multiplication, to be specliazed for chain types.
 */
template <typename VariVal, typename Vari1, typename Vari2, typename = void>
class multiply_vari {};

/**
 * Specialization of var multiplication for two `var_value`.
 */
template <typename VariVal, typename Vari1, typename Vari2>
class multiply_vari<VariVal, Vari1, Vari2, require_all_vari_t<Vari1, Vari2>>
    final : public op_vari<VariVal, Vari1*, Vari2*> {
  using op_vari<VariVal, Vari1*, Vari2*>::avi;
  using op_vari<VariVal, Vari1*, Vari2*>::bvi;

 public:
  multiply_vari(Vari1* avi, Vari2* bvi)
      : op_vari<VariVal, Vari1*, Vari2*>(avi->val_ * bvi->val_, avi, bvi) {}
  /**
   * `chain_impl` is called from `chain()` and exists so one specialized struct
   * can call either the scalar or matrix `chain()` methods. SFINAE only works
   * on "deduced" template types. So the trick here is to make template types,
   * T1 and T2, set their defaults to the class template types, then do the
   * regular `requires`. Since `chain_impl` has no inputs to deduce the
   * template types will always fall back to their default values. Since the
   * compiler has "deduced" these types we can these use the standard requires
   * to SFINAE out either the arithmetic or matrix version.
   */
  template <typename T1 = Vari1, typename T2 = Vari2,
            require_all_vari_vt<std::is_arithmetic, T1, T2>* = nullptr>
  inline void chain_impl() {
      avi()->adj_ += bvi()->val_ * this->adj_;
      bvi()->adj_ += avi()->val_ * this->adj_;
  }

  template <typename T1 = Vari1, typename T2 = Vari2,
            require_all_vari_vt<is_eigen, T1, T2>* = nullptr>
  inline void chain_impl() {
      avi()->adj_ += this->adj_ * bvi()->val_.transpose();
      bvi()->adj_ += avi()->val_.transpose() * this->adj_;
  }

  void chain() {
    if (unlikely(is_any_nan(avi()->val_, bvi()->val_))) {
      fill(avi()->adj_, NOT_A_NUMBER);
      fill(bvi()->adj_, NOT_A_NUMBER);
    } else {
     chain_impl();
    }
  }
};

/**
 * Specialization of var multiplication for `var_value` and arithmetic
 */
template <typename VariVal, typename Vari, typename Arith>
class multiply_vari<VariVal, Vari, Arith, require_vt_arithmetic<Arith>> final
    : public op_vari<VariVal, Vari*, Arith> {
  using op_vari<VariVal, Vari*, Arith>::avi;
  using op_vari<VariVal, Vari*, Arith>::bd;

 public:
  multiply_vari(Vari* avi, const Arith& b)
      : op_vari<VariVal, Vari*, Arith>(avi->val_ * b, avi, b) {}

  template <typename T1 = Vari, typename T2 = Arith,
            require_vari_vt<std::is_arithmetic, T1>* = nullptr,
            require_vt_arithmetic<T2>* = nullptr>
  inline void chain_impl() {
      avi()->adj_ += this->adj_ * bd();
  }

  template <typename T1 = Vari, typename T2 = Arith,
            require_vari_vt<is_eigen, T1>* = nullptr,
            require_vt_arithmetic<T2>* = nullptr>
  inline void chain_impl() {
      avi()->adj_ += this->adj_ * bd().transpose();
  }

  void chain() {
    if (unlikely(is_any_nan(avi()->val_, bd()))) {
      fill(avi()->adj_, NOT_A_NUMBER);
    } else {
      chain_impl();
    }
  }
};

/**
 * Specialization of var multiplication for arithmetic and `var_value`
 */
template <typename VariVal, typename Arith, typename Vari>
class multiply_vari<VariVal, Arith, Vari, require_vt_arithmetic<Arith>> final
    : public op_vari<VariVal, Arith, Vari*> {
  using op_vari<VariVal, Arith, Vari*>::ad;
  using op_vari<VariVal, Arith, Vari*>::bvi;

 public:
  multiply_vari(const Arith& a, Vari* bvi)
      : op_vari<VariVal, Arith, Vari*>(a * bvi->val_, a, bvi) {}

  template <typename T1 = Arith, typename T2 = Vari,
            require_vari_vt<std::is_arithmetic, T2>* = nullptr,
            require_vt_arithmetic<T1>* = nullptr>
  void chain_impl() {
    bvi()->adj_ += this->adj_ * ad();
  }

  template <typename T1 = Arith, typename T2 = Vari,
            require_vari_vt<is_eigen, T2>* = nullptr,
            require_vt_arithmetic<T1>* = nullptr>
  void chain_impl() {
    bvi()->adj_ += (this->adj_ * ad()).transpose();
  }

  void chain() {
    auto a = ad();
    auto b = bvi()->val_;
    if (unlikely(is_any_nan(b, a))) {
      fill(bvi()->adj_, NOT_A_NUMBER);
    } else {
      chain_impl();
    }
  }
};

/**
 * Deduces the return type for matrix multiplication of two types
 */
template <typename T1, typename T2, typename = void, typename = void>
struct mat_mul_return_type {};

// arithmetic is just double
template <typename T1, typename T2>
struct mat_mul_return_type<T1, T2, require_all_arithmetic_t<T1, T2>> {
  /**
   * FIXME: Should probs do something to promote to highest type given
   * something like float/double/int
   */
  using type = double;
};

struct mult_invoker {
  template <typename T1, typename T2>
  auto operator()(T1&& x, T2&& y) {
    return (x * y).eval();
  }
};

template <typename T1, typename T2>
struct mat_mul_return_type<T1, T2, require_any_eigen_t<T1, T2>> {
  using type = std::result_of_t<mult_invoker(T1, T2)>;
};
// helper alias
template <typename T1, typename T2>
using mat_mul_return_type_t = typename mat_mul_return_type<T1, T2>::type;
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
template <typename T1, typename T2, require_all_var_value_t<T1, T2>* = nullptr>
inline auto operator*(const T1& a, const T2& b) {
  using vari1 = get_var_vari_value_t<T1>;
  using vari2 = get_var_vari_value_t<T2>;
  using mat_return = internal::mat_mul_return_type_t<get_var_scalar_t<T1>, get_var_scalar_t<T2>>;
  using multiply_type = internal::multiply_vari<mat_return, vari1, vari2>;
  return var_value<mat_return>{new multiply_type(a.vi_, b.vi_)};
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
template <typename T, typename Arith, require_var_value_t<T>* = nullptr,
 require_vt_arithmetic<Arith>* = nullptr>
inline auto operator*(
    const T& a, const Arith& b) {
  using vari_type = get_var_vari_value_t<T>;
  using mat_return = internal::mat_mul_return_type_t<get_var_scalar_t<T>, Arith>;
  using multiply_type = internal::multiply_vari<mat_return, vari_type, Arith>;
  return var_value<mat_return>{new multiply_type(a.vi_, b)};
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
template <typename T, typename Arith, require_var_value_t<T>* = nullptr,
  require_vt_arithmetic<Arith>* = nullptr>
inline auto operator*(const Arith& a, const T& b) {
  using vari_type = get_var_vari_value_t<T>;
  using mat_return = internal::mat_mul_return_type_t<Arith, get_var_scalar_t<T>>;
  using multiply_type = internal::multiply_vari<mat_return, Arith, vari_type>;
  return var_value<mat_return>{new multiply_type(a, b.vi_)};
}

}  // namespace math
}  // namespace stan
#endif
