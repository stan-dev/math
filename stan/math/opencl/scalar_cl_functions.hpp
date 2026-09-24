#ifndef STAN_MATH_OPENCL_SCALAR_CL_FUNCTIONS_HPP
#define STAN_MATH_OPENCL_SCALAR_CL_FUNCTIONS_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/scalar_cl.hpp>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {
namespace opencl {

/** \ingroup opencl
 * Operators and element-wise functions whose operands are all scalars and at
 * least one is a device scalar. They are evaluated immediately on the device
 * with a single thread and return a device scalar. Expressions mixing device
 * scalars with matrices are built lazily by the kernel generator instead.
 *
 * These live only in `stan::math::opencl` and are found by argument dependent
 * lookup or by qualifying with `opencl::`.
 */
namespace internal {
/**
 * Checks if a type is a host arithmetic value or a device scalar holding one.
 */
template <typename T>
struct is_scalar_cl_operand
    : math::disjunction<std::is_arithmetic<std::decay_t<T>>,
                        is_prim_scalar_cl<T>> {};

/**
 * Enables a template if all types are host arithmetic values or device
 * scalars and at least one is a device scalar.
 */
template <typename... Types>
using require_scalar_cl_operands_t = require_t<
    math::conjunction<is_scalar_cl_operand<Types>...,
                      math::disjunction<is_prim_scalar_cl<Types>...>>>;
}  // namespace internal

#define STAN_OPENCL_SCALAR_CL_BINARY_OPERATOR(op)                          \
  template <typename T_a, typename T_b,                                    \
            internal::require_scalar_cl_operands_t<T_a, T_b>* = nullptr>   \
  inline ScalarCl<double> operator op(T_a && a, T_b && b) {                \
    return ScalarCl<double>(as_operation_cl(std::forward<T_a>(a))          \
                                op as_operation_cl(std::forward<T_b>(b))); \
  }

STAN_OPENCL_SCALAR_CL_BINARY_OPERATOR(+)
STAN_OPENCL_SCALAR_CL_BINARY_OPERATOR(-)
#undef STAN_OPENCL_SCALAR_CL_BINARY_OPERATOR

/**
 * Multiplies scalars, at least one of which is a device scalar.
 * @tparam T_a type of the first operand
 * @tparam T_b type of the second operand
 * @param a first operand
 * @param b second operand
 * @return device scalar holding the product
 */
template <typename T_a, typename T_b,
          internal::require_scalar_cl_operands_t<T_a, T_b>* = nullptr>
inline ScalarCl<double> operator*(T_a&& a, T_b&& b) {
  return ScalarCl<double>(
      stan::math::elt_multiply(as_operation_cl(std::forward<T_a>(a)),
                               as_operation_cl(std::forward<T_b>(b))));
}

/**
 * Divides scalars, at least one of which is a device scalar.
 * @tparam T_a type of the first operand
 * @tparam T_b type of the second operand
 * @param a first operand
 * @param b second operand
 * @return device scalar holding the quotient
 */
template <typename T_a, typename T_b,
          internal::require_scalar_cl_operands_t<T_a, T_b>* = nullptr>
inline ScalarCl<double> operator/(T_a&& a, T_b&& b) {
  return ScalarCl<double>(
      stan::math::elt_divide(as_operation_cl(std::forward<T_a>(a)),
                             as_operation_cl(std::forward<T_b>(b))));
}

/**
 * Negates a device scalar.
 * @tparam T type of the device scalar
 * @param a device scalar
 * @return device scalar holding the negation
 */
template <typename T, require_prim_scalar_cl_t<T>* = nullptr>
inline ScalarCl<double> operator-(T&& a) {
  return ScalarCl<double>(-as_operation_cl(std::forward<T>(a)));
}

#define STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(fun)               \
  template <typename T, require_prim_scalar_cl_t<T>* = nullptr> \
  inline ScalarCl<double> fun(T&& a) {                          \
    return ScalarCl<double>(                                    \
        stan::math::fun(as_operation_cl(std::forward<T>(a))));  \
  }

#define STAN_OPENCL_SCALAR_CL_BINARY_FUNCTION(fun)                       \
  template <typename T_a, typename T_b,                                  \
            internal::require_scalar_cl_operands_t<T_a, T_b>* = nullptr> \
  inline ScalarCl<double> fun(T_a&& a, T_b&& b) {                        \
    return ScalarCl<double>(                                             \
        stan::math::fun(as_operation_cl(std::forward<T_a>(a)),           \
                        as_operation_cl(std::forward<T_b>(b))));         \
  }

STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(rsqrt)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(sqrt)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(cbrt)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(exp)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(exp2)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(expm1)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(log)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(log2)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(log10)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(log1p)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(sin)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(sinh)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(cos)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(cosh)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(tan)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(tanh)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(asin)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(asinh)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(acos)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(acosh)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(atan)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(atanh)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(tgamma)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(lgamma)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(erf)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(erfc)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(floor)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(round)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(ceil)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(fabs)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(trunc)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(digamma)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(erfcx)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(log1m)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(log_inv_logit)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(log1m_exp)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(log1p_exp)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(inv_square)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(inv_logit)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(logit)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(Phi)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(Phi_approx)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(inv_Phi)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(log1m_inv_logit)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(trigamma)
STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION(square)

STAN_OPENCL_SCALAR_CL_BINARY_FUNCTION(fdim)
STAN_OPENCL_SCALAR_CL_BINARY_FUNCTION(fmax)
STAN_OPENCL_SCALAR_CL_BINARY_FUNCTION(fmin)
STAN_OPENCL_SCALAR_CL_BINARY_FUNCTION(fmod)
STAN_OPENCL_SCALAR_CL_BINARY_FUNCTION(hypot)
STAN_OPENCL_SCALAR_CL_BINARY_FUNCTION(ldexp)
STAN_OPENCL_SCALAR_CL_BINARY_FUNCTION(pow)
STAN_OPENCL_SCALAR_CL_BINARY_FUNCTION(copysign)
STAN_OPENCL_SCALAR_CL_BINARY_FUNCTION(beta)
STAN_OPENCL_SCALAR_CL_BINARY_FUNCTION(binomial_coefficient_log)
STAN_OPENCL_SCALAR_CL_BINARY_FUNCTION(lbeta)
STAN_OPENCL_SCALAR_CL_BINARY_FUNCTION(log_inv_logit_diff)
STAN_OPENCL_SCALAR_CL_BINARY_FUNCTION(log_diff_exp)
STAN_OPENCL_SCALAR_CL_BINARY_FUNCTION(multiply_log)
STAN_OPENCL_SCALAR_CL_BINARY_FUNCTION(lmultiply)

#undef STAN_OPENCL_SCALAR_CL_UNARY_FUNCTION
#undef STAN_OPENCL_SCALAR_CL_BINARY_FUNCTION

}  // namespace opencl
}  // namespace math
}  // namespace stan

#endif
#endif
