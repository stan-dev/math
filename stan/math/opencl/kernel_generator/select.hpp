#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_SELECT_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_SELECT_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/is_valid_expression.hpp>
#include <stan/math/opencl/kernel_generator/common_return_scalar.hpp>
#include <set>
#include <string>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {

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
  using base::var_name;

 protected:
  using base::arguments_;

 public:
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
   * generates kernel code for this (select) operation.
   * @param i row index variable name
   * @param j column index variable name
   * @param var_name_condition variable name of the condition expression
   * @param var_name_else variable name of the then expression
   * @param var_name_then variable name of the else expression
   * @return part of kernel with code for this expression
   */
  inline kernel_parts generate(const std::string& i, const std::string& j,
                               const std::string& var_name_condition,
                               const std::string& var_name_then,
                               const std::string& var_name_else) const {
    kernel_parts res{};
    res.body = type_str<Scalar>() + " " + var_name + " = " + var_name_condition
               + " ? " + var_name_then + " : " + var_name_else + ";\n";
    return res;
  }

  /**
   * View of a matrix that would be the result of evaluating this expression.
   * @return view
   */
  inline matrix_cl_view view() const {
    matrix_cl_view condition_view = std::get<0>(arguments_).view();
    matrix_cl_view then_view = std::get<1>(arguments_).view();
    matrix_cl_view else_view = std::get<2>(arguments_).view();
    return both(either(then_view, else_view), both(condition_view, then_view));
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
template <typename T_condition, typename T_then, typename T_else,
          typename
          = require_all_valid_expressions_t<T_condition, T_then, T_else>>
inline select_<as_operation_cl_t<T_condition>, as_operation_cl_t<T_then>,
               as_operation_cl_t<T_else>>
select(T_condition&& condition, T_then&& then, T_else&& els) {  // NOLINT
  return {as_operation_cl(std::forward<T_condition>(condition)),
          as_operation_cl(std::forward<T_then>(then)),
          as_operation_cl(std::forward<T_else>(els))};
}

}  // namespace math
}  // namespace stan
#endif
#endif
