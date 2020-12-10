#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_ELT_FUNCTION_CL_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_ELT_FUNCTION_CL_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/kernels/device_functions/digamma.hpp>
#include <stan/math/opencl/kernels/device_functions/log1m_exp.hpp>
#include <stan/math/opencl/kernels/device_functions/log1m_inv_logit.hpp>
#include <stan/math/opencl/kernels/device_functions/log1p_exp.hpp>
#include <stan/math/opencl/kernels/device_functions/logit.hpp>
#include <stan/math/opencl/kernels/device_functions/inv_logit.hpp>
#include <stan/math/opencl/kernels/device_functions/inv_square.hpp>
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
    static_assert(std::is_floating_point<                                      \
                      typename std::remove_reference_t<T>::Scalar>::value,     \
                  #fun                                                         \
                  ": all arguments must be expression with floating point "    \
                  "return type!");                                             \
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
ADD_UNARY_FUNCTION_WITH_INCLUDES(log1m_exp,
                                 opencl_kernels::log1m_exp_device_function)
ADD_UNARY_FUNCTION_WITH_INCLUDES(log1p_exp,
                                 opencl_kernels::log1p_exp_device_function)
ADD_UNARY_FUNCTION_WITH_INCLUDES(inv_square,
                                 opencl_kernels::inv_square_device_function)
ADD_UNARY_FUNCTION_WITH_INCLUDES(inv_logit,
                                 opencl_kernels::inv_logit_device_function)
ADD_UNARY_FUNCTION_WITH_INCLUDES(logit, opencl_kernels::logit_device_function)
ADD_UNARY_FUNCTION_WITH_INCLUDES(
    log1m_inv_logit, opencl_kernels::log1m_inv_logit_device_function)

ADD_CLASSIFICATION_FUNCTION(isfinite, {-rows() + 1, cols() - 1})
ADD_CLASSIFICATION_FUNCTION(isinf,
                            this->template get_arg<0>().extreme_diagonals())
ADD_CLASSIFICATION_FUNCTION(isnan,
                            this->template get_arg<0>().extreme_diagonals())

#undef ADD_UNARY_FUNCTION_WITH_INCLUDES
#undef ADD_UNARY_FUNCTION
#undef ADD_UNARY_FUNCTION_PASS_ZERO
#undef ADD_CLASSIFICATION_FUNCTION

/** @}*/
}  // namespace math
}  // namespace stan
#endif
#endif
