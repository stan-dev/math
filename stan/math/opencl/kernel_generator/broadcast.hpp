#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_BROADCAST_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_BROADCAST_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/is_valid_expression.hpp>
#include <string>
#include <type_traits>
#include <set>
#include <utility>

namespace stan {
namespace math {

/**
 * Represents a broadcasting operation_cl in kernel generator expressions.
 * @tparam T type of argument
 * @tparam Rows whether to broadcast rows
 * @tparam Cols whether to broadcast columns
 */
template <typename T, bool Rows, bool Cols>
class broadcast_
    : public operation_cl<broadcast_<T, Rows, Cols>,
                          typename std::remove_reference_t<T>::Scalar, T> {
 public:
  using Scalar = typename std::remove_reference_t<T>::Scalar;
  using base = operation_cl_cl<broadcast_<T, Rows, Cols>, Scalar, T>;
  using base::instance;
  using base::var_name;

  /**
   * Constructor
   * @param a expression
   */
  explicit broadcast_(T&& a) : base(std::forward<T>(a)) {
    const char* function = "broadcast";
    if (Rows) {
      check_size_match(function, "Rows of ", "a", a.rows(), "", "", 1);
    }
    if (Cols) {
      check_size_match(function, "Columns of ", "a", a.cols(), "", "", 1);
    }
  }

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() {
    auto&& arg_copy = std::get<0>(arguments_).deep_copy();
    return block_<std::remove_reference_t<decltype(arg_copy)>>{
        std::move(arg_copy), start_row_, start_col_, rows_, cols_};
  }

  /**
   * Generates kernel code for this and nested expressions.
   * @param[in,out] generated set of already generated operations
   * @param ng name generator for this kernel
   * @param i row index variable name
   * @param j column index variable name
   * @return part of kernel with code for this and nested expressions
   */
  inline kernel_parts generate(const std::string& i, const std::string& j,
                               const std::string& var_name_arg) const {
    kernel_parts res
        = a_.generate(generated, ng, Rows ? "0" : i, Cols ? "0" : j,
                      Rows ? "1" : rows, Cols ? "1" : rows);
    var_name = a_.var_name;
    return res;
  }

  /**
   * Sets kernel arguments for this and nested expressions.
   * @param[in,out] generated set of expressions that already set their kernel
   * arguments
   * @param kernel kernel to set arguments on
   * @param[in,out] arg_num consecutive number of the first argument to set.
   * This is incremented for each argument set by this function.
   */
  inline void set_args(std::set<int>& generated, cl::Kernel& kernel,
                       int& arg_num) const {
    a_.set_args(generated, kernel, arg_num);
  }

  /**
   * Adds event for any matrices used by this or nested expressions.
   * @param e the event to add
   */
  inline void add_event(cl::Event& e) const { a_.add_event(e); }

  /**
   * Number of rows of a matrix that would be the result of evaluating this
   * expression.
   * @return number of rows
   */
  inline int rows() const { return Rows ? base::dynamic : a_.rows(); }

  /**
   * Number of columns of a matrix that would be the result of evaluating this
   * expression.
   * @return number of columns
   */
  inline int cols() const { return Cols ? base::dynamic : a_.cols(); }

  /**
   * View of a matrix that would be the result of evaluating this expression.
   * @return view
   */
  inline matrix_cl_view view() const {
    matrix_cl_view view = a_.view();
    if (Rows) {
      view = either(view, matrix_cl_view::Lower);
    }
    if (Cols) {
      view = either(view, matrix_cl_view::Upper);
    }
    return view;
  }
};

/**
 * Brodcast an expression in specified dimension(s). If broadcasting a
 * dimension, that dimension of the input must be equal to 1. Further
 * expressions can use this expression as if had any size in broadcast
 * dimension, repeating the values.
 * @tparam Rows whether to broadcast rows
 * @tparam Cols whether to broadcast rows
 * @tparam T type of input expression
 * @param a input expression
 * @return broadcast expression
 */
template <bool Rows, bool Cols, typename T,
          typename = enable_if_all_valid_expressions_and_none_scalar<T>>
inline broadcast_<as_operation_cl_t<T>, Rows, Cols> broadcast(T&& a) {
  return broadcast_<as_operation_cl_t<T>, Rows, Cols>(
      as_operation_cl(std::forward<T>(a)));
}

}  // namespace math
}  // namespace stan
#endif
#endif
