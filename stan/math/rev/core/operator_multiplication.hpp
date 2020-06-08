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

template <typename T>
struct is_not_var : bool_constant<!is_var<T>::value> {};

namespace internal {

/**
 * Deduces the return type for matrix multiplication of two types
 */
template <typename T1, typename T2, typename = void>
struct mat_mul_return_type {};

// helper alias
template <typename T1, typename T2>
using mat_mul_return_t = std::decay_t<
    typename mat_mul_return_type<std::decay_t<T1>, std::decay_t<T2>>::type>;

// arithmetic is just double
template <typename T1, typename T2>
struct mat_mul_return_type<T1, T2, require_all_arithmetic_t<T1, T2>> {
  using type = double;
};

template <typename T1, typename T2>
struct mat_mul_return_type<T1, T2, require_any_eigen_t<T1, T2>> {
  using type = decltype((std::declval<T1>() * std::declval<T2>()).eval());
};

template <typename T1, typename T2, typename = void>
struct mat_mul_var_return_type {};

template <typename T1, typename T2>
struct mat_mul_var_return_type<T1, T2, require_all_var_t<T1, T2>> {
  using type = mat_mul_return_t<value_type_t<T1>, value_type_t<T2>>;
};

template <typename T1, typename T2>
struct mat_mul_var_return_type<T1, T2,
                               require_all_t<is_var<T1>, is_not_var<T2>>> {
  using type = mat_mul_return_t<value_type_t<T1>, T2>;
};

template <typename T1, typename T2>
struct mat_mul_var_return_type<T1, T2,
                               require_all_t<is_not_var<T1>, is_var<T2>>> {
  using type = mat_mul_return_t<T1, value_type_t<T2>>;
};

template <typename T1, typename T2>
using mat_mul_var_return_t =
    typename mat_mul_var_return_type<std::decay_t<T1>, std::decay_t<T2>>::type;

/**
 * Base for multiplication, to be specliazed for chain types.
 */
template <typename Var1, typename Var2, typename = void>
class multiply_vari;

/**
 * Specialization of var multiplication for two `var_value`.
 */
template <typename Var1, typename Var2>
class multiply_vari<Var1, Var2, require_all_var_t<Var1, Var2>> final
    : public op_vari<mat_mul_var_return_t<Var1, Var2>, Var1, Var2> {
  using op_vari_mul = op_vari<mat_mul_var_return_t<Var1, Var2>, Var1, Var2>;
  using op_vari_mul::avi;
  using op_vari_mul::bvi;

 public:
  using op_vari_mul::return_t;
  template <typename T1, typename T2>
  multiply_vari(T1* avi, T2* bvi)
      : op_vari<mat_mul_var_return_t<Var1, Var2>, Var1, Var2>(
            avi->val_ * bvi->val_, avi, bvi) {}
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
  template <typename T1 = Var1, typename T2 = Var2,
            require_all_var_vt<std::is_arithmetic, T1, T2>* = nullptr>
  inline void chain_impl() {
    avi()->adj_ += bvi()->val_ * this->adj_;
    bvi()->adj_ += avi()->val_ * this->adj_;
  }

  template <typename T1 = Var1, typename T2 = Var2,
            require_all_var_vt<is_eigen, T1, T2>* = nullptr>
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
template <typename Var, typename Arith>
class multiply_vari<Var, Arith, require_all_t<is_var<Var>, is_not_var<Arith>>>
    final : public op_vari<mat_mul_var_return_t<Var, Arith>, Var, Arith> {
  using op_vari_mul = op_vari<mat_mul_var_return_t<Var, Arith>, Var, Arith>;
  using op_vari_mul::avi;
  using op_vari_mul::bd;

 public:
  using op_vari_mul::return_t;
  template <typename T1, typename T2>
  multiply_vari(T1* avi, const T2& b) : op_vari_mul(avi->val_ * b, avi, b) {}

  template <typename T1 = Var, typename T2 = Arith,
            require_var_vt<std::is_arithmetic, T1>* = nullptr,
            require_arithmetic_t<T2>* = nullptr>
  inline void chain_impl() {
    avi()->adj_ += this->adj_ * bd();
  }

  template <typename T1 = Var, typename T2 = Arith,
            require_var_vt<is_eigen, T1>* = nullptr,
            require_vt_arithmetic<T2>* = nullptr>
  inline void chain_impl() {
    avi()->adj_ += this->adj_ * bd().transpose();
  }
  // NOTE: THIS IS WRONG
  template <typename T1 = Arith, typename T2 = Var,
            require_eigen_t<T1>* = nullptr,
            require_var_vt<std::is_arithmetic, T2>* = nullptr>
  void chain_impl() {
    avi()->adj_ += this->adj_.sum();
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
template <typename Arith, typename Var>
class multiply_vari<Arith, Var, require_all_t<is_not_var<Arith>, is_var<Var>>>
    final : public op_vari<mat_mul_var_return_t<Arith, Var>, Arith, Var> {
  using op_vari_mul = op_vari<mat_mul_var_return_t<Arith, Var>, Arith, Var>;
  using op_vari_mul::ad;
  using op_vari_mul::bvi;

 public:
  using op_vari_mul::return_t;
  template <typename T2>
  multiply_vari(const Arith& a, T2* bvi) : op_vari_mul(a * bvi->val_, a, bvi) {}

  template <typename T1 = Arith, typename T2 = Var,
            require_arithmetic_t<T1>* = nullptr,
            require_var_vt<std::is_arithmetic, T2>* = nullptr>
  void chain_impl() {
    bvi()->adj_ += this->adj_ * ad();
  }

  template <typename T1 = Arith, typename T2 = Var,
            require_vt_arithmetic<T1>* = nullptr,
            require_var_vt<is_eigen, T2>* = nullptr>
  void chain_impl() {
    bvi()->adj_ += (this->adj_ * ad()).transpose();
  }
  // NOTE: THIS IS WRONG
  template <typename T1 = Arith, typename T2 = Var,
            require_eigen_t<T1>* = nullptr,
            require_var_vt<std::is_arithmetic, T2>* = nullptr>
  void chain_impl() {
    bvi()->adj_ += this->adj_.sum();
  }

  void chain() {
    if (unlikely(is_any_nan(bvi()->val_, ad()))) {
      fill(bvi()->adj_, NOT_A_NUMBER);
    } else {
      chain_impl();
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
template <typename T1, typename T2, require_all_var_value_t<T1, T2>* = nullptr>
inline auto operator*(const T1& a, const T2& b) {
  using multiply_type = internal::multiply_vari<T1, T2>;
  using mat_return = typename multiply_type::return_t;
  return var_value<mat_return>(new multiply_type(a.vi_, b.vi_));
}

/**
 * idk how to name this, but we need something that at SFINAE
 * that allows mixes of stan scalar types or var<eig> with eigen but not
 * var and Eig<double> which should differ to the old implimentation for dynamic
 * types.
 */
template <typename T, typename S, typename = void>
struct is_conformable : std::true_type {};

template <typename T, typename S>
struct is_conformable<T, S, require_var_t<T>>
    : bool_constant<!is_eigen<S>::value> {};

template <typename T, typename S>
using require_conformable_t = require_t<is_conformable<T, S>>;

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
          require_st_arithmetic<Arith>* = nullptr>
//,
//          require_conformable_t<T, Arith>* = nullptr>
inline auto operator*(const T& a, const Arith& b) {
  using multiply_type = internal::multiply_vari<T, Arith>;
  using mat_return = typename multiply_type::return_t;
  return var_value<mat_return>(new multiply_type(a.vi_, b));
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
          require_st_arithmetic<Arith>* = nullptr>
//,
//          require_conformable_t<T, Arith>* = nullptr>
inline auto operator*(const Arith& a, const T& b) {
  using multiply_type = internal::multiply_vari<Arith, T>;
  using mat_return = typename multiply_type::return_t;
  return var_value<mat_return>(new multiply_type(a, b.vi_));
}

}  // namespace math
}  // namespace stan
#endif
