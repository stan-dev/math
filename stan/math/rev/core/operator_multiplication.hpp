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
template <typename VariVal, typename Vari1, typename Vari2, typename = void>
class multiply_vari {};

template <typename VariVal, typename Vari1, typename Vari2>
class multiply_vari<VariVal, Vari1, Vari2, require_all_vari_t<Vari1, Vari2>>
    final : public op_vari<VariVal, Vari1*, Vari2*> {
  using op_vari<VariVal, Vari1*, Vari2*>::avi;
  using op_vari<VariVal, Vari1*, Vari2*>::bvi;

 public:
  multiply_vari(Vari1* avi, Vari2* bvi)
      : op_vari<VariVal, Vari1*, Vari2*>(avi->val_ * bvi->val_, avi, bvi) {}
  template <typename T1 = Vari1, typename T2 = Vari2,
   require_all_vari_vt<std::is_arithmetic, T1, T2>* = nullptr>
  void chain_impl() {
    if (unlikely(is_any_nan(avi()->val_, bvi()->val_))) {
      avi()->adj_ = NOT_A_NUMBER;
      bvi()->adj_ = NOT_A_NUMBER;
    } else {
      avi()->adj_ += bvi()->val_ * this->adj_;
      bvi()->adj_ += avi()->val_ * this->adj_;
    }
  }

  template <typename T1 = Vari1, typename T2 = Vari2,
   require_all_vari_vt<is_eigen, T1, T2>* = nullptr>
  void chain_impl() {
    if (unlikely(is_any_nan(avi()->val_, bvi()->val_))) {
      avi()->adj_.fill(NOT_A_NUMBER);
      bvi()->adj_.fill(NOT_A_NUMBER);
    } else {
      avi()->adj_ += this->adj_ * bvi()->val_.transpose();
      bvi()->adj_ += avi()->val_.transpose() * this->adj_;
    }
  }

  void chain() {
    chain_impl();
  }
};

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
   require_arithmetic_t<T2>* = nullptr>
  void chain_impl() {
    if (unlikely(is_any_nan(avi()->val_, bd()))) {
      avi()->adj_ = NOT_A_NUMBER;
    } else {
      avi()->adj_ += this->adj_ * bd();
    }
  }

  template <typename T1 = Vari, typename T2 = Arith,
   require_vari_vt<is_eigen, T1>* = nullptr,
   require_vt_arithmetic<T2>* = nullptr>
  void chain_impl() {
    if (unlikely(is_any_nan(avi()->val_, bd()))) {
      avi()->adj_.fill(NOT_A_NUMBER);
    } else {
      avi()->adj_ += this->adj_ * bd().transpose();
    }
  }

  void chain() {
    chain_impl();
  }
};

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
   require_arithmetic_t<T1>* = nullptr>
  void chain_impl() {
    if (unlikely(is_any_nan(bvi()->val_, ad()))) {
      bvi()->adj_ = NOT_A_NUMBER;
    } else {
      bvi()->adj_ += this->adj_ * ad();
    }
  }

  template <typename T1 = Arith, typename T2 = Vari,
   require_vari_vt<is_eigen, T2>* = nullptr,
   require_vt_arithmetic<T1>* = nullptr>
  void chain_impl() {
    if (unlikely(is_any_nan(bvi()->val_, ad()))) {
      bvi()->adj_.fill(NOT_A_NUMBER);
    } else {
      bvi()->adj_ += (this->adj_ * ad()).transpose();
    }
  }

  void chain() {
    chain_impl();
  }
};

template <typename T1, typename T2, typename = void, typename = void>
struct mat_mul_return_type {};

template <typename T1, typename T2>
struct mat_mul_return_type<T1, T2, require_all_arithmetic_t<T1, T2>> {
  using type = double;
};

template <typename T1, typename T2>
struct mat_mul_return_type<T1, T2, require_all_eigen_matrix_t<T1, T2>> {
  using type = std::decay_t<T1>;
};

template <typename T1, typename T2>
struct mat_mul_return_type<T1, T2, require_eigen_matrix_t<T1>, require_eigen_col_vector_t<T2>> {
  using type = std::decay_t<T2>;
};

template <typename T1, typename T2>
struct mat_mul_return_type<T1, T2, require_eigen_row_vector_t<T1>, require_eigen_matrix_t<T2>> {
  using type = std::decay_t<T1>;
};

template <typename T1, typename T2>
struct mat_mul_return_type<T1, T2, require_eigen_row_vector_t<T1>, require_eigen_col_vector_t<T2>> {
  using type = value_type_t<T1>;
};

template <typename T1, typename T2>
struct mat_mul_return_type<T1, T2, require_eigen_col_vector_t<T1>, require_eigen_row_vector_t<T2>> {
  using type = Eigen::Matrix<value_type_t<T1>, -1, -1>;
};

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
template <typename T1, typename T2>
inline var_value<internal::mat_mul_return_type_t<T1, T2>> operator*(const var_value<T1>& a, const var_value<T2>& b) {
  using mat_return = internal::mat_mul_return_type_t<T1, T2>;
  return {new internal::multiply_vari<mat_return,
     vari_value<T1>, vari_value<T2>>(a.vi_, b.vi_)};
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
inline var_value<internal::mat_mul_return_type_t<T, Arith>> operator*(const var_value<T>& a, const Arith& b) {
  if (is_all_equal(b, 1.0)) {
    return a;
  }
  using mat_return = internal::mat_mul_return_type_t<T, Arith>;
  return {new internal::multiply_vari<mat_return, vari_value<T>, Arith>(a.vi_, b)};
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
inline var_value<internal::mat_mul_return_type_t<Arith, T>> operator*(const Arith& a, const var_value<T>& b) {
  if (is_all_equal(a, 1.0)) {
    return b;
  }
  using mat_return = internal::mat_mul_return_type_t<Arith, T>;
  return {new internal::multiply_vari<mat_return, Arith, vari_value<T>>(a, b.vi_)};
}

}  // namespace math
}  // namespace stan
#endif
