#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_ELT_FUNCTION_CL_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_ELT_FUNCTION_CL_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/kernels/device_functions/binomial_coefficient_log.hpp>
#include <stan/math/opencl/kernels/device_functions/beta.hpp>
#include <stan/math/opencl/kernels/device_functions/digamma.hpp>
#include <stan/math/opencl/kernels/device_functions/inv_logit.hpp>
#include <stan/math/opencl/kernels/device_functions/inv_Phi.hpp>
#include <stan/math/opencl/kernels/device_functions/inv_square.hpp>
#include <stan/math/opencl/kernels/device_functions/lbeta.hpp>
#include <stan/math/opencl/kernels/device_functions/lgamma_stirling.hpp>
#include <stan/math/opencl/kernels/device_functions/lgamma_stirling_diff.hpp>
#include <stan/math/opencl/kernels/device_functions/lmultiply.hpp>
#include <stan/math/opencl/kernels/device_functions/log_inv_logit.hpp>
#include <stan/math/opencl/kernels/device_functions/log_inv_logit_diff.hpp>
#include <stan/math/opencl/kernels/device_functions/log_diff_exp.hpp>
#include <stan/math/opencl/kernels/device_functions/log1m.hpp>
#include <stan/math/opencl/kernels/device_functions/log1m_exp.hpp>
#include <stan/math/opencl/kernels/device_functions/log1m_inv_logit.hpp>
#include <stan/math/opencl/kernels/device_functions/log1p_exp.hpp>
#include <stan/math/opencl/kernels/device_functions/logit.hpp>
#include <stan/math/opencl/kernels/device_functions/multiply_log.hpp>
#include <stan/math/opencl/kernels/device_functions/Phi.hpp>
#include <stan/math/opencl/kernels/device_functions/Phi_approx.hpp>
#include <stan/math/opencl/kernels/device_functions/trigamma.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernel_generator/common_return_scalar.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <array>
#include <string>
#include <type_traits>
#include <set>
#include <utility>

namespace stan {
namespace math {

/** \addtogroup opencl_kernel_generator
 *  @{
 */

/**
 * Represents an element-wise function in kernel generator expressions.
 * @tparam Derived derived type
 * @tparam T type of argument
 * @tparam Scal type of the scalar of result
 */
template <typename Derived, typename Scal, typename... T>
class elt_function_cl : public operation_cl<Derived, Scal, T...> {
 public:
  using Scalar = Scal;
  using base = operation_cl<Derived, Scalar, T...>;
  using base::var_name_;

  /**
   * Constructor
   * @param fun function
   * @param args argument expression(s)
   */
  elt_function_cl(const std::string& fun, T&&... args)  // NOLINT
      : base(std::forward<T>(args)...), fun_(fun) {}

  /**
   * Generates kernel code for this expression.
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @param view_handled whether whether caller already handled matrix view
   * @param var_names_arg variable names of the nested expressions
   * @return part of kernel with code for this expression
   */
  inline kernel_parts generate(
      const std::string& row_index_name, const std::string& col_index_name,
      const bool view_handled,
      std::conditional_t<false, T, const std::string&>... var_names_arg) const {
    kernel_parts res{};

    for (const char* incl : base::derived().includes) {
      res.includes += incl;
    }
    std::array<std::string, sizeof...(T)> var_names_arg_arr
        = {(var_names_arg + ", ")...};
    std::string var_names_list = std::accumulate(
        var_names_arg_arr.begin(), var_names_arg_arr.end(), std::string());
    res.body = type_str<Scalar>() + " " + var_name_ + " = " + fun_ + "((double)"
               + var_names_list.substr(0, var_names_list.size() - 2) + ");\n";
    return res;
  }

 protected:
  std::string fun_;
};

/**
 * Generates a class and function for a general binary function that is defined
 * by OpenCL or in the included code.
 * @param fun function
 * @param ... function sources to include into kernel
 */
#define ADD_BINARY_FUNCTION_WITH_INCLUDES(fun, ...)                         \
  template <typename T1, typename T2>                                       \
  class fun##_ : public elt_function_cl<fun##_<T1, T2>, double, T1, T2> {   \
    using base = elt_function_cl<fun##_<T1, T2>, double, T1, T2>;           \
    using base::arguments_;                                                 \
                                                                            \
   public:                                                                  \
    using base::rows;                                                       \
    using base::cols;                                                       \
    static const std::vector<const char*> includes;                         \
    explicit fun##_(T1&& a, T2&& b)                                         \
        : base(#fun, std::forward<T1>(a), std::forward<T2>(b)) {            \
      if (a.rows() != base::dynamic && b.rows() != base::dynamic) {         \
        check_size_match(#fun, "Rows of ", "a", a.rows(), "rows of ", "b",  \
                         b.rows());                                         \
      }                                                                     \
      if (a.cols() != base::dynamic && b.cols() != base::dynamic) {         \
        check_size_match(#fun, "Columns of ", "a", a.cols(), "columns of ", \
                         "b", b.cols());                                    \
      }                                                                     \
    }                                                                       \
    inline auto deep_copy() const {                                         \
      auto&& arg1_copy = this->template get_arg<0>().deep_copy();           \
      auto&& arg2_copy = this->template get_arg<1>().deep_copy();           \
      return fun##_<std::remove_reference_t<decltype(arg1_copy)>,           \
                    std::remove_reference_t<decltype(arg2_copy)>>{          \
          std::move(arg1_copy), std::move(arg2_copy)};                      \
    }                                                                       \
    inline std::pair<int, int> extreme_diagonals() const {                  \
      return {-rows() + 1, cols() - 1};                                     \
    }                                                                       \
  };                                                                        \
                                                                            \
  template <typename T1, typename T2,                                       \
            require_all_kernel_expressions_t<T1, T2>* = nullptr,            \
            require_any_not_stan_scalar_t<T1, T2>* = nullptr>               \
  inline fun##_<as_operation_cl_t<T1>, as_operation_cl_t<T2>> fun(T1&& a,   \
                                                                  T2&& b) { \
    return fun##_<as_operation_cl_t<T1>, as_operation_cl_t<T2>>(            \
        as_operation_cl(std::forward<T1>(a)),                               \
        as_operation_cl(std::forward<T2>(b)));                              \
  }                                                                         \
  template <typename T1, typename T2>                                       \
  const std::vector<const char*> fun##_<T1, T2>::includes{__VA_ARGS__};

/**
 * Generates a class and function for a general unary function that is defined
 * by OpenCL or in the included code.
 * @param fun function
 * @param ... function sources to include into kernel
 */
#define ADD_UNARY_FUNCTION_WITH_INCLUDES(fun, ...)                             \
  template <typename T>                                                        \
  class fun##_ : public elt_function_cl<fun##_<T>, double, T> {                \
    using base = elt_function_cl<fun##_<T>, double, T>;                        \
    using base::arguments_;                                                    \
                                                                               \
   public:                                                                     \
    using base::rows;                                                          \
    using base::cols;                                                          \
    static const std::vector<const char*> includes;                            \
    explicit fun##_(T&& a) : base(#fun, std::forward<T>(a)) {}                 \
    inline auto deep_copy() const {                                            \
      auto&& arg_copy = this->template get_arg<0>().deep_copy();               \
      return fun##_<std::remove_reference_t<decltype(arg_copy)>>{              \
          std::move(arg_copy)};                                                \
    }                                                                          \
    inline std::pair<int, int> extreme_diagonals() const {                     \
      return {-rows() + 1, cols() - 1};                                        \
    }                                                                          \
  };                                                                           \
                                                                               \
  template <typename T, typename Cond                                          \
                        = require_all_kernel_expressions_and_none_scalar_t<T>> \
  inline fun##_<as_operation_cl_t<T>> fun(T&& a) {                             \
    return fun##_<as_operation_cl_t<T>>(as_operation_cl(std::forward<T>(a)));  \
  }                                                                            \
  template <typename T>                                                        \
  const std::vector<const char*> fun##_<T>::includes{__VA_ARGS__};

/**
 * Generates a class and function for a general unary function that is defined
 * by OpenCL.
 * @param fun function
 */
#define ADD_UNARY_FUNCTION(fun) ADD_UNARY_FUNCTION_WITH_INCLUDES(fun)

/**
 * Generates a class and function for an unary function, defined by OpenCL with
 * special property that it passes trough zero. That is \f$ f(0)=0 \f$. Such a
 * function can have triangular view equal to its argument's.
 * @param fun function name
 */
#define ADD_UNARY_FUNCTION_PASS_ZERO(fun)                                      \
  template <typename T>                                                        \
  class fun##_ : public elt_function_cl<fun##_<T>, double, T> {                \
    using base = elt_function_cl<fun##_<T>, double, T>;                        \
    using base::arguments_;                                                    \
                                                                               \
   public:                                                                     \
    using base::rows;                                                          \
    using base::cols;                                                          \
    static constexpr auto view_transitivness = std::make_tuple(true);          \
    static const std::vector<const char*> includes;                            \
    explicit fun##_(T&& a) : base(#fun, std::forward<T>(a)) {}                 \
    inline auto deep_copy() const {                                            \
      auto&& arg_copy = this->template get_arg<0>().deep_copy();               \
      return fun##_<std::remove_reference_t<decltype(arg_copy)>>{              \
          std::move(arg_copy)};                                                \
    }                                                                          \
  };                                                                           \
                                                                               \
  template <typename T, typename Cond                                          \
                        = require_all_kernel_expressions_and_none_scalar_t<T>> \
  inline fun##_<as_operation_cl_t<T>> fun(T&& a) {                             \
    return fun##_<as_operation_cl_t<T>>(as_operation_cl(std::forward<T>(a)));  \
  }                                                                            \
  template <typename T>                                                        \
  const std::vector<const char*> fun##_<T>::includes{};

/**
 * Generates a class and function for a classification function, defined by
 * OpenCL.
 * @param fun function name
 * @param ... code for determining extreme diagonals
 */
#define ADD_CLASSIFICATION_FUNCTION(fun, ...)                                  \
  template <typename T>                                                        \
  class fun##_ : public elt_function_cl<fun##_<T>, bool, T> {                  \
    using base = elt_function_cl<fun##_<T>, bool, T>;                          \
    using base::arguments_;                                                    \
                                                                               \
   public:                                                                     \
    using base::rows;                                                          \
    using base::cols;                                                          \
    static constexpr auto view_transitivness = std::make_tuple(true);          \
    static const std::vector<const char*> includes;                            \
    explicit fun##_(T&& a) : base(#fun, std::forward<T>(a)) {}                 \
    inline auto deep_copy() const {                                            \
      auto&& arg_copy = this->template get_arg<0>().deep_copy();               \
      return fun##_<std::remove_reference_t<decltype(arg_copy)>>{              \
          std::move(arg_copy)};                                                \
    }                                                                          \
    inline std::pair<int, int> extreme_diagonals() const {                     \
      return __VA_ARGS__;                                                      \
    }                                                                          \
  };                                                                           \
                                                                               \
  template <typename T, typename Cond                                          \
                        = require_all_kernel_expressions_and_none_scalar_t<T>> \
  inline fun##_<as_operation_cl_t<T>> fun(T&& a) {                             \
    return fun##_<as_operation_cl_t<T>>(as_operation_cl(std::forward<T>(a)));  \
  }                                                                            \
  template <typename T>                                                        \
  const std::vector<const char*> fun##_<T>::includes{};

ADD_UNARY_FUNCTION(rsqrt)
ADD_UNARY_FUNCTION_PASS_ZERO(sqrt)
ADD_UNARY_FUNCTION_PASS_ZERO(cbrt)

ADD_UNARY_FUNCTION(exp)
ADD_UNARY_FUNCTION(exp2)
ADD_UNARY_FUNCTION_PASS_ZERO(expm1)

ADD_UNARY_FUNCTION(log)
ADD_UNARY_FUNCTION(log2)
ADD_UNARY_FUNCTION(log10)
ADD_UNARY_FUNCTION_PASS_ZERO(log1p)

ADD_UNARY_FUNCTION_PASS_ZERO(sin)
ADD_UNARY_FUNCTION_PASS_ZERO(sinh)
ADD_UNARY_FUNCTION(cos)
ADD_UNARY_FUNCTION(cosh)
ADD_UNARY_FUNCTION_PASS_ZERO(tan)
ADD_UNARY_FUNCTION_PASS_ZERO(tanh)
ADD_UNARY_FUNCTION_PASS_ZERO(asin)
ADD_UNARY_FUNCTION_PASS_ZERO(asinh)
ADD_UNARY_FUNCTION(acos)
ADD_UNARY_FUNCTION(acosh)
ADD_UNARY_FUNCTION_PASS_ZERO(atan)
ADD_UNARY_FUNCTION_PASS_ZERO(atanh)

ADD_UNARY_FUNCTION(tgamma)
ADD_UNARY_FUNCTION(lgamma)
ADD_UNARY_FUNCTION_PASS_ZERO(erf)
ADD_UNARY_FUNCTION(erfc)

ADD_UNARY_FUNCTION_PASS_ZERO(floor)
ADD_UNARY_FUNCTION_PASS_ZERO(round)
ADD_UNARY_FUNCTION_PASS_ZERO(ceil)
ADD_UNARY_FUNCTION_PASS_ZERO(fabs)
ADD_UNARY_FUNCTION_PASS_ZERO(trunc)

ADD_UNARY_FUNCTION_WITH_INCLUDES(digamma,
                                 opencl_kernels::digamma_device_function)
ADD_UNARY_FUNCTION_WITH_INCLUDES(log1m, opencl_kernels::log1m_device_function)
ADD_UNARY_FUNCTION_WITH_INCLUDES(log_inv_logit,
                                 opencl_kernels::log_inv_logit_device_function)
ADD_UNARY_FUNCTION_WITH_INCLUDES(log1m_exp,
                                 opencl_kernels::log1m_exp_device_function)
ADD_UNARY_FUNCTION_WITH_INCLUDES(log1p_exp,
                                 opencl_kernels::log1p_exp_device_function)
ADD_UNARY_FUNCTION_WITH_INCLUDES(inv_square,
                                 opencl_kernels::inv_square_device_function)
ADD_UNARY_FUNCTION_WITH_INCLUDES(inv_logit,
                                 opencl_kernels::inv_logit_device_function)
ADD_UNARY_FUNCTION_WITH_INCLUDES(logit, opencl_kernels::logit_device_function)
ADD_UNARY_FUNCTION_WITH_INCLUDES(Phi, opencl_kernels::phi_device_function)
ADD_UNARY_FUNCTION_WITH_INCLUDES(Phi_approx,
                                 opencl_kernels::inv_logit_device_function,
                                 opencl_kernels::phi_approx_device_function)
ADD_UNARY_FUNCTION_WITH_INCLUDES(inv_Phi, opencl_kernels::log1m_device_function,
                                 opencl_kernels::phi_device_function,
                                 opencl_kernels::inv_phi_device_function)
ADD_UNARY_FUNCTION_WITH_INCLUDES(
    log1m_inv_logit, opencl_kernels::log1m_inv_logit_device_function)
ADD_UNARY_FUNCTION_WITH_INCLUDES(trigamma,
                                 opencl_kernels::trigamma_device_function)
ADD_UNARY_FUNCTION_WITH_INCLUDES(
    square,
    "\n#ifndef STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_SQUARE\n"
    "#define STAN_MATH_OPENCL_KERNELS_DEVICE_FUNCTIONS_SQUARE\n"
    "double square(double x){return x*x;}\n"
    "#endif\n")

ADD_CLASSIFICATION_FUNCTION(isfinite, {-rows() + 1, cols() - 1})
ADD_CLASSIFICATION_FUNCTION(isinf,
                            this->template get_arg<0>().extreme_diagonals())
ADD_CLASSIFICATION_FUNCTION(isnan,
                            this->template get_arg<0>().extreme_diagonals())

ADD_BINARY_FUNCTION_WITH_INCLUDES(fdim)
ADD_BINARY_FUNCTION_WITH_INCLUDES(fmax)
ADD_BINARY_FUNCTION_WITH_INCLUDES(fmin)
ADD_BINARY_FUNCTION_WITH_INCLUDES(fmod)
ADD_BINARY_FUNCTION_WITH_INCLUDES(hypot)
ADD_BINARY_FUNCTION_WITH_INCLUDES(ldexp)
ADD_BINARY_FUNCTION_WITH_INCLUDES(pow)
ADD_BINARY_FUNCTION_WITH_INCLUDES(copysign)

ADD_BINARY_FUNCTION_WITH_INCLUDES(
    beta, stan::math::opencl_kernels::beta_device_function)
ADD_BINARY_FUNCTION_WITH_INCLUDES(
    binomial_coefficient_log,
    stan::math::opencl_kernels::lgamma_stirling_device_function,
    stan::math::opencl_kernels::lgamma_stirling_diff_device_function,
    stan::math::opencl_kernels::lbeta_device_function,
    stan::math::opencl_kernels::binomial_coefficient_log_device_function)
ADD_BINARY_FUNCTION_WITH_INCLUDES(
    lbeta, stan::math::opencl_kernels::lgamma_stirling_device_function,
    stan::math::opencl_kernels::lgamma_stirling_diff_device_function,
    stan::math::opencl_kernels::lbeta_device_function)
ADD_BINARY_FUNCTION_WITH_INCLUDES(
    log_inv_logit_diff, opencl_kernels::log1p_exp_device_function,
    opencl_kernels::log1m_exp_device_function,
    opencl_kernels::log_inv_logit_diff_device_function)
ADD_BINARY_FUNCTION_WITH_INCLUDES(log_diff_exp,
                                  opencl_kernels::log1m_exp_device_function,
                                  opencl_kernels::log_diff_exp_device_function)
ADD_BINARY_FUNCTION_WITH_INCLUDES(
    multiply_log, stan::math::opencl_kernels::multiply_log_device_function)
ADD_BINARY_FUNCTION_WITH_INCLUDES(
    lmultiply, stan::math::opencl_kernels::lmultiply_device_function)

#undef ADD_BINARY_FUNCTION_WITH_INCLUDES
#undef ADD_UNARY_FUNCTION_WITH_INCLUDES
#undef ADD_UNARY_FUNCTION
#undef ADD_UNARY_FUNCTION_PASS_ZERO
#undef ADD_CLASSIFICATION_FUNCTION

/** @}*/
}  // namespace math
}  // namespace stan
#endif
#endif
