#ifndef STAN_MATH_OPENCL_REV_MULTIPLY_HPP
#define STAN_MATH_OPENCL_REV_MULTIPLY_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator/is_kernel_expression.hpp>
#include <stan/math/opencl/multiply.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>
#include <stan/math/prim/fun/value_of.hpp>

namespace stan {
namespace math {

namespace internal {
  template <typename T>
  struct matrix_cl_alloc : public chainable_alloc {
    T mat_;
    matrix_cl_alloc(const T& x) : mat_(x) {}
  };
  template <typename T1, typename T2, typename T3, typename = void, typename = void>
  struct multiply_cl_vari;

  template <typename T1, typename T2, typename T3>
  struct multiply_cl_vari<T1, T2, T3, require_var_t<T1>, require_not_var_t<T2>> : public vari {
    T1 a_;
    matrix_cl_alloc<T2>* b_;
    T3 res_;
    multiply_cl_vari(const T1& a, const T2& b, const T3& c) : vari(0),
      a_(a), b_(new matrix_cl_alloc<T2>(b)), res_(c) {}

    void chain() {
      puts("Got 1");
      a_.adj() = a_.adj() + res_.adj() * transpose(b_->mat_);
    }

  };

  template <typename T1, typename T2, typename T3>
  struct multiply_cl_vari<T1, T2, T3, require_not_var_t<T1>, require_var_t<T2>> : public vari {
    matrix_cl_alloc<T1>* a_;
    T2 b_;
    T3 res_;
    multiply_cl_vari(const T1& a, const T2& b, const T3& c) : vari(0),
      a_(new matrix_cl_alloc<std::decay_t<T1>>(a)), b_(b), res_(c) {}

    void chain() {
      puts("Got 2");
      b_.adj() = b_.adj() + transpose(a_->mat_) * res_.adj();
    }

  };

}

/**
 * Matrix multiplication of two reverse mode matrices and/or kernel generator
 * expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 * @return Matrix product of given arguments
 */
template <
    typename T_a, typename T_b,
    require_all_nonscalar_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
    require_all_var_t<T_a, T_b>* = nullptr>
inline auto multiply(const T_a& a, const T_b& b) {
  check_size_match("multiply ((OpenCL))", "A.cols()", a.cols(), "B.rows()",
                   b.rows());

  var_value<matrix_cl<double>> res = value_of(b) * value_of(a);

  reverse_pass_callback([a, b, res]() mutable {
      a.adj() = a.adj() + res.adj() * transpose(b.val());
      b.adj() = b.adj() + transpose(a.val()) * res.adj();
  });
  return res;
}

template <
    typename T_a, typename T_b,
    require_all_nonscalar_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
    require_not_var_t<T_a>* = nullptr, require_var_t<T_b>* = nullptr>
inline auto multiply(const T_a& a, const T_b& b) {
  check_size_match("multiply ((OpenCL))", "A.cols()", a.cols(), "B.rows()",
                   b.rows());

  var_value<matrix_cl<double>> res = a * b.val();
  auto* do_mult = new internal::multiply_cl_vari<plain_type_t<T_a>,
   plain_type_t<T_b>, var_value<matrix_cl<double>>>(a, b, res);
  return res;
}

template <
    typename T_a, typename T_b,
    require_all_nonscalar_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
    require_var_t<T_a>* = nullptr, require_not_var_t<T_b>* = nullptr>
inline auto multiply(const T_a& a, const T_b& b) {
  check_size_match("multiply ((OpenCL))", "A.cols()", a.cols(), "B.rows()",
                   b.rows());
  var_value<matrix_cl<double>> res = a.val() * b;
  auto* do_mult = new internal::multiply_cl_vari<plain_type_t<T_a>,
  plain_type_t<T_b>, var_value<matrix_cl<double>>>(a, b, res);
  return res;
}

/**
 * Matrix multiplication of two reverse mode matrices and/or kernel generator
 * expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 * @return Matrix product of given arguments
 */
template <
    typename T_a, typename T_b,
    require_all_nonscalar_prim_or_rev_kernel_expression_t<T_a, T_b>* = nullptr,
    require_any_var_t<T_a, T_b>* = nullptr>
inline auto operator*(const T_a& a, const T_b& b) {
  return multiply(a, b);
}

}  // namespace math
}  // namespace stan

#endif
#endif
