#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_ROWWISE_REDUCTION_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_ROWWISE_REDUCTION_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/is_valid_expression.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <set>
#include <string>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {
/**
 * Represents a rowwise reduction in kernel generator expressions.
 * @tparam Derived derived type
 * @tparam T type of first argument
 * @tparam operation type with member function generate that accepts two
 * variable names and returns OpenCL source code for reduction operation_cl
 * @tparam PassZero whether \c operation passes trough zeros
 * @tparam Rowwise whether this is row wise reduction
 * @tparam Colwise whether this is column wise reduction
 */

template <typename Derived, typename T, typename operation, bool PassZero>
class rowwise_reduction
    : public operation_cl<Derived, typename std::remove_reference_t<T>::Scalar,
                          T> {
 public:
  using Scalar = typename std::remove_reference_t<T>::Scalar;
  using base = operation_cl<Derived, Scalar, T>;
  using base::var_name;

 protected:
  std::string init_;
  using base::arguments_;

 public:
  /**
   * Constructor
   * @param a the expression to reduce
   * @param init OpenCL source code of initialization value for reduction
   */
  rowwise_reduction(T&& a, const std::string& init)
      : base(std::forward<T>(a)), init_(init) {}

  /**
   * Generates kernel code for this expression.
   * @param i row index variable name
   * @param j column index variable name
   * @param var_name_arg name of the variable in kernel that holds argument to
   * this expression
   * @return part of kernel with code for this expression
   */
  inline kernel_parts generate(const std::string& i, const std::string& j,
                               const std::string& var_name_arg) const {
    kernel_parts res;
    res.body_prefix
        = type_str<Scalar>() + " " + var_name + " = " + init_ + ";\n";
    if (PassZero) {
      res.body_prefix += "for(int " + var_name + "_j = contains_nonzero("
                         + var_name + "_view, LOWER) ? 0 : " + i + "; "
                         + var_name + "_j < (contains_nonzero(" + var_name
                         + "_view, UPPER) ? " + var_name + "_cols : min("
                         + var_name + "_cols, " + i + " + 1)); " + var_name
                         + "_j++){\n";
    } else {
      res.body_prefix += "for(int " + var_name + "_j = 0; " + var_name + "_j < "
                         + var_name + "_cols; " + var_name + "_j++){\n";
    }
    res.body += var_name + " = " + operation::generate(var_name, var_name_arg)
                + ";\n}\n";
    res.args = "int " + var_name + "_view, int " + var_name + "_cols, ";
    return res;
  }

  /**
   * Sets offset of block to indices of the argument expression
   * @param[in, out] i row index
   * @param[in, out] j column index
   */
  inline void modify_argument_indices(std::string& i, std::string& j) const {
    j = var_name + "_j";
  }

  /**
   * Sets kernel arguments for this and nested expressions.
   * @param[in,out] generated set of expressions that already set their kernel
   * arguments
   * @param kernel kernel to set arguments on
   * @param[in,out] arg_num consecutive number of the first argument to set.
   * This is incremented for each argument set by this function.
   */
  inline void set_args(std::set<const operation_cl_base*>& generated,
                       cl::Kernel& kernel, int& arg_num) const {
    if (generated.count(this) == 0) {
      generated.insert(this);
      std::get<0>(arguments_).set_args(generated, kernel, arg_num);
      kernel.setArg(arg_num++, std::get<0>(arguments_).view());
      kernel.setArg(arg_num++, std::get<0>(arguments_).cols());
    }
  }

  /**
   * Number of columns of a matrix that would be the result of evaluating this
   * expression.
   * @return number of columns
   */
  inline int cols() const { return 1; }

  /**
   * View of a matrix that would be the result of evaluating this expression.
   * @return view
   */
  matrix_cl_view view() const { return matrix_cl_view::Entire; }

  /**
   * Determine index of top diagonal written.
   * @return number of columns
   */
  inline int top_diagonal() const { return 1; }
};

/**
 * Operation for sum reduction.
 */
struct sum_op {
  /**
   * Generates sum reduction kernel code.
   * @param a first variable
   * @param b second variable
   * @return reduction code
   */
  inline static std::string generate(const std::string& a,
                                     const std::string& b) {
    return a + " + " + b;
  }
};

/**
 * Represents rowwise sum reduction in kernel generator expressions.
 * @tparam T type of expression
 */
template <typename T>
class rowwise_sum_
    : public rowwise_reduction<rowwise_sum_<T>, T, sum_op, true> {
 public:
  explicit rowwise_sum_(T&& a)
      : rowwise_reduction<rowwise_sum_<T>, T, sum_op, true>(std::forward<T>(a),
                                                            "0") {}
};

/**
 * Rowwise sum reduction of a kernel generator expression.
 * @tparam T type of input expression
 * @param a expression to reduce
 * @return sum
 */
template <typename T,
          typename = require_all_valid_expressions_and_none_scalar_t<T>>
inline rowwise_sum_<as_operation_cl_t<T>> rowwise_sum(T&& a) {
  return rowwise_sum_<as_operation_cl_t<T>>(
      as_operation_cl(std::forward<T>(a)));
}

/**
 * Operation for max reduction.
 * @tparam T type to reduce
 */
template <typename T>
struct max_op {
  /**
   * Generates max reduction kernel code.
   * @param a first variable
   * @param b second variable
   * @return reduction code
   */
  inline static std::string generate(const std::string& a,
                                     const std::string& b) {
    if (std::is_floating_point<T>()) {
      return "fmax(" + a + ", " + b + ")";
    }
    return "max(" + a + ", " + b + ")";
  }

  inline static std::string init() {
    if (std::is_floating_point<T>()) {
      return "-INFINITY";
    }
    return "INT_MIN";
  }
};

/**
 * Represents rowwise max reduction in kernel generator expressions.
 * @tparam T type of expression
 */
template <typename T>
class rowwise_max_
    : public rowwise_reduction<
          rowwise_max_<T>, T,
          max_op<typename std::remove_reference_t<T>::Scalar>, false> {
 public:
  using op = max_op<typename std::remove_reference_t<T>::Scalar>;
  explicit rowwise_max_(T&& a)
      : rowwise_reduction<rowwise_max_<T>, T, op, false>(std::forward<T>(a),
                                                         op::init()) {}
};

/**
 * Rowwise max reduction of a kernel generator expression.
 * @tparam T type of input expression
 * @param a expression to reduce
 * @return max
 */
template <typename T,
          typename = require_all_valid_expressions_and_none_scalar_t<T>>
inline rowwise_max_<as_operation_cl_t<T>> rowwise_max(T&& a) {
  return rowwise_max_<as_operation_cl_t<T>>(
      as_operation_cl(std::forward<T>(a)));
}

/**
 * Operation for min reduction.
 * @tparam T type to reduce
 */
template <typename T>
struct min_op {
  /**
   * Generates min reduction kernel code.
   * @param a first variable
   * @param b second variable
   * @return reduction code
   */
  inline static std::string generate(const std::string& a,
                                     const std::string& b) {
    if (std::is_floating_point<T>()) {
      return "fmin(" + a + ", " + b + ")";
    }
    return "min(" + a + ", " + b + ")";
  }

  inline static std::string init() {
    if (std::is_floating_point<T>()) {
      return "INFINITY";
    }
    return "INT_MAX";
  }
};

/**
 * Represents rowwise min reduction in kernel generator expressions.
 * @tparam T type of expression
 */
template <typename T>
class rowwise_min_
    : public rowwise_reduction<
          rowwise_min_<T>, T,
          min_op<typename std::remove_reference_t<T>::Scalar>, false> {
 public:
  using op = min_op<typename std::remove_reference_t<T>::Scalar>;
  explicit rowwise_min_(T&& a)
      : rowwise_reduction<rowwise_min_<T>, T, op, false>(std::forward<T>(a),
                                                         op::init()) {}
};

/**
 * Min reduction of a kernel generator expression.
 * @tparam T type of input expression
 * @param a expression to reduce
 * @return min
 */
template <typename T,
          typename = require_all_valid_expressions_and_none_scalar_t<T>>
inline rowwise_min_<as_operation_cl_t<T>> rowwise_min(T&& a) {
  return rowwise_min_<as_operation_cl_t<T>>(
      as_operation_cl(std::forward<T>(a)));
}

}  // namespace math
}  // namespace stan

#endif
#endif
