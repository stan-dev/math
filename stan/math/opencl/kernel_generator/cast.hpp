#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_CAST_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_CAST_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
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
 * Represents a typecast os scalar in kernel generator expressions.
 * @tparam Derived derived type
 * @tparam T type of argument
 * @tparam Scal type of the scalar of result
 */
template <typename Scal, typename T>
class cast_ : public operation_cl<cast_<Scal, T>, Scal, T> {
 public:
  using Scalar = Scal;
  using base = operation_cl<cast_<Scal, T>, Scalar, T>;
  using base::var_name_;

  /**
   * Constructor
   * @param args argument expression(s)
   */
  explicit cast_(T&& arg) : base(std::forward<T>(arg)) {}

  /**
   * Generates kernel code for this expression.
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @param view_handled whether whether caller already handled matrix view
   * @param var_names_arg variable names of the nested expressions
   * @return part of kernel with code for this expression
   */
  inline kernel_parts generate(const std::string& row_index_name,
                               const std::string& col_index_name,
                               const bool view_handled,
                               const std::string& var_name_arg) const {
    kernel_parts res{};

    res.body = type_str<Scalar>() + " " + var_name_ + " = ("
               + type_str<Scalar>() + ")" + var_name_arg + ";\n";
    return res;
  }

  inline auto deep_copy() const {
    auto&& arg_copy = this->template get_arg<0>().deep_copy();
    return cast_<Scalar, std::remove_reference_t<decltype(arg_copy)>>{
        std::move(arg_copy)};
  }
};

/**
 * Typecast a kernel generator expression scalar.
 *
 * @tparam T type of argument
 * @param a input argument
 * @return Typecast of given expression
 */
template <typename Scalar, typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline auto cast(T&& a) {
  auto&& a_operation = as_operation_cl(std::forward<T>(a)).deep_copy();
  return cast_<Scalar, std::remove_reference_t<decltype(a_operation)>>(
      std::move(a_operation));
}

/**
 * Typecast a scalar.
 *
 * @tparam T type of argument
 * @param a input argument
 * @return Typecast of given expression
 */
template <typename Scalar, typename T, require_stan_scalar_t<T>* = nullptr>
inline Scalar cast(T a) {
  return a;
}

/** @}*/
}  // namespace math
}  // namespace stan
#endif
#endif
