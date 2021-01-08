#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_OPTIONALBROADCAST_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_OPTIONALBROADCAST_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <limits>
#include <string>
#include <type_traits>
#include <map>
#include <utility>

namespace stan {
namespace math {

/**
 * Represents an optional broadcasting operation in kernel generator
 * expressions.
 * @tparam T type of argument
 * @tparam Colwise whether to broadcast colwise
 * @tparam Rowwise whether to broadcast rowwise
 */
template <typename T, bool Colwise, bool Rowwise>
class optional_broadcast_
    : public operation_cl<optional_broadcast_<T, Colwise, Rowwise>,
                          typename std::remove_reference_t<T>::Scalar, T> {
 public:
  using Scalar = typename std::remove_reference_t<T>::Scalar;
  using base
      = operation_cl<optional_broadcast_<T, Colwise, Rowwise>, Scalar, T>;
  using base::var_name_;

  /**
   * Constructor
   * @param a expression
   */
  explicit optional_broadcast_(T&& a) : base(std::forward<T>(a)) {}

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& arg_copy = this->template get_arg<0>().deep_copy();
    return optional_broadcast_<std::remove_reference_t<decltype(arg_copy)>,
                               Colwise, Rowwise>{std::move(arg_copy)};
  }

  /**
   * Generates kernel code for this and nested expressions.
   * @param row_idx_name  row index variable name
   * @param col_idx_name  column index variable name
   * @param view_handled whether whether caller already handled matrix view
   * @param var_name_arg name of the variable in kernel that holds argument to
   * this expression
   * @return part of kernel with code for this and nested expressions
   */
  inline kernel_parts generate(const std::string& row_idx_name,
                               const std::string& col_idx_name,
                               const bool view_handled,
                               const std::string& var_name_arg) const {
    kernel_parts res;
    res.body
        += type_str<Scalar>() + " " + var_name_ + " = " + var_name_arg + ";\n";
    if (Colwise) {
      res.args += "int " + var_name_ + "is_multirow, ";
    }
    if (Rowwise) {
      res.args += "int " + var_name_ + "is_multicol, ";
    }
    return res;
  }

  /**
   * Sets index/indices along broadcasted dimmension(s) to 0.
   * @param[in, out] row_idx_name row index
   * @param[in, out] col_idx_name  column index
   */
  inline void modify_argument_indices(std::string& row_idx_name,
                                      std::string& col_idx_name) const {
    if (Colwise) {
      row_idx_name = "(" + row_idx_name + " * " + var_name_ + "is_multirow)";
    }
    if (Rowwise) {
      col_idx_name = "(" + col_idx_name + " * " + var_name_ + "is_multicol)";
    }
  }

  /**
   * Sets kernel arguments for this and nested expressions.
   * @param[in,out] generated map from (pointer to) already generated local
   * operations to variable names
   * @param[in,out] generated_all map from (pointer to) already generated all
   * operations to variable names
   * @param kernel kernel to set arguments on
   * @param[in,out] arg_num consecutive number of the first argument to set.
   * This is incremented for each argument set by this function.
   */
  inline void set_args(std::map<const void*, const char*>& generated,
                       std::map<const void*, const char*>& generated_all,
                       cl::Kernel& kernel, int& arg_num) const {
    if (generated.count(this) == 0) {
      generated[this] = "";
      std::map<const void*, const char*> generated2;
      this->template get_arg<0>().set_args(generated2, generated_all, kernel,
                                           arg_num);
      if (Colwise) {
        kernel.setArg(arg_num++, static_cast<int>(
                                     this->template get_arg<0>().rows() != 1));
      }
      if (Rowwise) {
        kernel.setArg(arg_num++, static_cast<int>(
                                     this->template get_arg<0>().cols() != 1));
      }
    }
  }

  /**
   * Number of rows of a matrix that would be the result of evaluating this
   * expression.
   * @return number of rows
   */
  inline int rows() const {
    return Colwise && this->template get_arg<0>().rows() == 1
               ? base::dynamic
               : this->template get_arg<0>().rows();
  }

  /**
   * Number of columns of a matrix that would be the result of evaluating this
   * expression.
   * @return number of columns
   */
  inline int cols() const {
    return Rowwise && this->template get_arg<0>().cols() == 1
               ? base::dynamic
               : this->template get_arg<0>().cols();
  }

  /**
   * View of a matrix that would be the result of evaluating this expression.
   * @return view
   */
  inline matrix_cl_view view() const {
    matrix_cl_view view = this->template get_arg<0>().view();
    if (Colwise && this->template get_arg<0>().rows() == 1) {
      view = either(view, matrix_cl_view::Lower);
    }
    if (Rowwise && this->template get_arg<0>().cols() == 1) {
      view = either(view, matrix_cl_view::Upper);
    }
    return view;
  }

  /**
   * Determine index of bottom diagonal written.
   * @return index of bottom diagonal
   */
  inline int bottom_diagonal() const {
    if (Colwise && this->template get_arg<0>().rows() == 1) {
      return std::numeric_limits<int>::min();
    } else {
      return this->template get_arg<0>().bottom_diagonal();
    }
  }

  /**
   * Determine index of top diagonal written.
   * @return index of top diagonal
   */
  inline int top_diagonal() const {
    if (Rowwise && this->template get_arg<0>().cols() == 1) {
      return std::numeric_limits<int>::max();
    } else {
      return this->template get_arg<0>().top_diagonal();
    }
  }
};

/**
 * Broadcast an expression in specified dimension(s) if the size along that
 * dimension equals 1. In that case further
 * expressions can use this expression as if it had any size in the broadcast
 * dimension, repeating the values.
 *
 * Broadcasting evaluates the expression argument multiple times. For
 * performance reasons don't broadcast slow operations. Instead evaluate them in
 * a separate kernel.
 * @tparam Colwise whether to broadcast colwise
 * @tparam Rowwise whether to broadcast rowwise
 * @tparam T type of input expression
 * @param a input expression
 * @return broadcast expression
 */
template <bool Colwise, bool Rowwise, typename T,
          typename = require_all_kernel_expressions_and_none_scalar_t<T>>
inline optional_broadcast_<as_operation_cl_t<T>, Colwise, Rowwise>
optional_broadcast(T&& a) {
  auto&& a_operation = as_operation_cl(std::forward<T>(a)).deep_copy();
  return optional_broadcast_<as_operation_cl_t<T>, Colwise, Rowwise>(
      std::move(a_operation));
}

/**
 * Broadcast an expression in rowwise dimmension if the number of columns equals
 * to 1. In that case further expressions can use this expression as if had any
 * number of columns, repeating the values.
 *
 * Broadcasting evaluates argument expression multiple times. For performance
 * reasons don't broadcast slow operations. Instead evaluate them in a separate
 * kernel.
 * @tparam T type of input expression
 * @param a input expression
 * @return broadcast expression
 */
template <typename T,
          typename = require_all_kernel_expressions_and_none_scalar_t<T>>
inline auto rowwise_optional_broadcast(T&& a) {
  return optional_broadcast<false, true>(std::forward<T>(a));
}

/**
 * Broadcast an expression in colwise dimmension if the number of rows equals
 * to 1. In that case further expressions can use this expression as if had any
 * number of rows, repeating the values.
 *
 * Broadcasting evaluates argument expression multiple times. For performance
 * reasons don't broadcast slow operations. Instead evaluate them in a separate
 * kernel.
 * @tparam T type of input expression
 * @param a input expression
 * @return broadcast expression
 */
template <typename T,
          typename = require_all_kernel_expressions_and_none_scalar_t<T>>
inline auto colwise_optional_broadcast(T&& a) {
  return optional_broadcast<true, false>(std::forward<T>(a));
}

}  // namespace math
}  // namespace stan

#endif
#endif
