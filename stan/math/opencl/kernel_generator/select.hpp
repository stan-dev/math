#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_SELECT_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_SELECT_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/common_return_scalar.hpp>
#include <algorithm>
#include <set>
#include <string>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {

/** \addtogroup opencl_kernel_generator
 *  @{
 */

/**
 * Represents a selection operation in kernel generator expressions. This is
 * element wise ternary operator <code>condition ? then : els</code>, also
 * equivalent to Eigen's \c .select().
 * @tparam Derived derived type
 * @tparam T_condition type of condition
 * @tparam T_then type of then expression
 * @tparam T_else type of else expression
 */
template <typename T_condition, typename T_then, typename T_else>
class select_ : public operation_cl<select_<T_condition, T_then, T_else>,
                                    common_scalar_t<T_then, T_else>,
                                    T_condition, T_then, T_else> {
 public:
  using Scalar = common_scalar_t<T_then, T_else>;
  using base = operation_cl<select_<T_condition, T_then, T_else>, Scalar,
                            T_condition, T_then, T_else>;
  using base::var_name_;

  /**
   * Constructor
   * @param condition condition expression
   * @param then then expression
   * @param els else expression
   */
  select_(T_condition&& condition, T_then&& then, T_else&& els)  // NOLINT
      : base(std::forward<T_condition>(condition), std::forward<T_then>(then),
             std::forward<T_else>(els)) {
    if (condition.rows() != base::dynamic && then.rows() != base::dynamic) {
      check_size_match("select", "Rows of ", "condition", condition.rows(),
                       "rows of ", "then", then.rows());
    }
    if (condition.cols() != base::dynamic && then.cols() != base::dynamic) {
      check_size_match("select", "Columns of ", "condition", condition.cols(),
                       "columns of ", "then", then.cols());
    }

    if (condition.rows() != base::dynamic && els.rows() != base::dynamic) {
      check_size_match("select", "Rows of ", "condition", condition.rows(),
                       "rows of ", "else", els.rows());
    }
    if (condition.cols() != base::dynamic && els.cols() != base::dynamic) {
      check_size_match("select", "Columns of ", "condition", condition.cols(),
                       "columns of ", "else", els.cols());
    }
  }

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& condition_copy = this->template get_arg<0>().deep_copy();
    auto&& then_copy = this->template get_arg<1>().deep_copy();
    auto&& else_copy = this->template get_arg<2>().deep_copy();
    return select_<std::remove_reference_t<decltype(condition_copy)>,
                   std::remove_reference_t<decltype(then_copy)>,
                   std::remove_reference_t<decltype(else_copy)>>(
        std::move(condition_copy), std::move(then_copy), std::move(else_copy));
  }

  /**
   * Generates kernel code for this (select) operation.
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @param view_handled whether whether caller already handled matrix view
   * @param var_name_condition variable name of the condition expression
   * @param var_name_else variable name of the then expression
   * @param var_name_then variable name of the else expression
   * @return part of kernel with code for this expression
   */
  inline kernel_parts generate(const std::string& row_index_name,
                               const std::string& col_index_name,
                               const bool view_handled,
                               const std::string& var_name_condition,
                               const std::string& var_name_then,
                               const std::string& var_name_else) const {
    kernel_parts res{};
    res.body = type_str<Scalar>() + " " + var_name_ + " = " + var_name_condition
               + " ? " + var_name_then + " : " + var_name_else + ";\n";
    return res;
  }

  /**
   * Determine indices of extreme sub- and superdiagonals written.
   * @return pair of indices - bottom and top diagonal
   */
  inline std::pair<int, int> extreme_diagonals() const {
    using std::max;
    using std::min;
    std::pair<int, int> condition_diags
        = this->template get_arg<0>().extreme_diagonals();
    std::pair<int, int> then_diags
        = this->template get_arg<1>().extreme_diagonals();
    std::pair<int, int> else_diags
        = this->template get_arg<2>().extreme_diagonals();
    // Where the condition is 0 we get else's values. Otherwise we get the more
    // extreme of then's and else's.
    return {max(min(then_diags.first, else_diags.first),
                min(condition_diags.first, else_diags.first)),
            min(max(then_diags.second, else_diags.second),
                max(condition_diags.second, else_diags.second))};
  }
};

/**
 * Selection operation on kernel generator expressions. This is element wise
 * ternary operator <code> condition ? then : els </code>.
 * @tparam T_condition type of condition expression
 * @tparam T_then type of then expression
 * @tparam T_else type of else expression
 * @param condition condition expression
 * @param then then expression
 * @param els else expression
 * @return selection operation expression
 */
template <
    typename T_condition, typename T_then, typename T_else,
    require_all_kernel_expressions_t<T_condition, T_then, T_else>* = nullptr,
    require_any_not_arithmetic_t<T_condition, T_then, T_else>* = nullptr>
inline select_<as_operation_cl_t<T_condition>, as_operation_cl_t<T_then>,
               as_operation_cl_t<T_else>>
select(T_condition&& condition, T_then&& then, T_else&& els) {  // NOLINT
  return {as_operation_cl(std::forward<T_condition>(condition)),
          as_operation_cl(std::forward<T_then>(then)),
          as_operation_cl(std::forward<T_else>(els))};
}

/**
 * Scalar overload of the selection operation.
 * @tparam T_then type of then scalar
 * @tparam T_else type of else scalar
 * @param condition condition
 * @param then then result
 * @param els else result
 * @return `condition ? then : els`
 */
template <typename T_then, typename T_else,
          require_all_arithmetic_t<T_then, T_else>* = nullptr>
inline std::common_type_t<T_then, T_else> select(bool condition, T_then then,
                                                 T_else els) {
  return condition ? then : els;
}

/** @}*/
}  // namespace math
}  // namespace stan
#endif
#endif
