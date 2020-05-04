#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_ROWWISE_REDUCTION_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_ROWWISE_REDUCTION_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/broadcast.hpp>
#include <stan/math/opencl/kernel_generator/binary_operation.hpp>
#include <stan/math/opencl/kernel_generator/is_valid_kernel_expression.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <set>
#include <string>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {
namespace internal {

/**
 * Implementation of an optimization for usage of rowwise reduction in
 * matrix-vector multiplication.
 */
template <typename Arg>
struct matvec_mul_opt {
  // in general the optimization is not possible
  enum { is_possible = 0 };

  static matrix_cl_view view(const Arg&) { return matrix_cl_view::Entire; }

  static kernel_parts get_kernel_parts(
      const Arg& a, std::set<const operation_cl_base*>& generated,
      name_generator& name_gen, const std::string& i, const std::string& j) {
    return {};
  }
};

template <typename Mat, typename VecT>
struct matvec_mul_opt<
    elewise_multiplication_<Mat, broadcast_<VecT, true, false>>> {
  // if the argument of rowwise reduction is multiplication with a broadcast
  // vector we can do the optimization
  enum { is_possible = 1 };
  using Arg = elewise_multiplication_<Mat, broadcast_<VecT, true, false>>;

  /**
   * Return view of the vector.
   * @param a argument to rowwise reduction (multiplication with second factor
   * being broadcast vector)
   * @return view
   */
  static matrix_cl_view view(const Arg& a) {
    return a.template get_arg<1>().template get_arg<0>().view();
  }

  /**
   * Generates kernel code for the argument of rowwise reduction, applying the
   * optimization - ignoring the triangular view of the vector, as it is already
   * handeled by rowwise reduction.
   * @param mul argument of the rowwise reduction
   * @param[in,out] generated set of (pointer to) already generated operations
   * @param name_gen name generator for this kernel
   * @param i row index variable name
   * @param j column index variable name
   * @return part of kernel with code for this and nested expressions
   */
  static kernel_parts get_kernel_parts(
      const Arg& mul, std::set<const operation_cl_base*>& generated,
      name_generator& name_gen, const std::string& i, const std::string& j) {
    kernel_parts res{};
    if (generated.count(&mul) == 0) {
      mul.var_name = name_gen.generate();
      generated.insert(&mul);

      const auto& matrix = mul.template get_arg<0>();
      const auto& broadcast = mul.template get_arg<1>();
      res = matrix.get_kernel_parts(generated, name_gen, i, j, true);
      if (generated.count(&broadcast) == 0) {
        broadcast.var_name = name_gen.generate();
        generated.insert(&broadcast);

        const auto& vec = broadcast.template get_arg<0>();
        std::string i_bc = i;
        std::string j_bc = j;
        broadcast.modify_argument_indices(i_bc, j_bc);
        res += vec.get_kernel_parts(generated, name_gen, i_bc, j_bc, true);
        res += broadcast.generate(i, j, true, vec.var_name);
      }
      res += mul.generate(i, j, true, matrix.var_name, broadcast.var_name);
    }
    return res;
  }
};

}  // namespace internal

/** \addtogroup opencl_kernel_generator
 *  @{
 */

/**
 * Represents a rowwise reduction in kernel generator expressions.
 * @tparam Derived derived type
 * @tparam T type of first argument
 * @tparam operation type with member function generate that accepts two
 * variable names and returns OpenCL source code for reduction operation_cl
 * @tparam PassZero whether \c operation passes trough zeros
 */
template <typename Derived, typename T, typename operation, bool PassZero>
class rowwise_reduction
    : public operation_cl<Derived, typename std::remove_reference_t<T>::Scalar,
                          T> {
 public:
  using T_no_ref = std::remove_reference_t<T>;
  using Scalar = typename T_no_ref::Scalar;
  using base = operation_cl<Derived, Scalar, T>;
  using base::var_name;

 protected:
  std::string init_;

 public:
  using base::rows;
  /**
   * Constructor
   * @param a the expression to reduce
   * @param init OpenCL source code of initialization value for reduction
   */
  explicit rowwise_reduction(T&& a, const std::string& init)
      : base(std::forward<T>(a)), init_(init) {}

  /**
   * Generates kernel code for this and nested expressions.
   * @param[in,out] generated set of (pointer to) already generated operations
   * @param name_gen name generator for this kernel
   * @param i row index variable name
   * @param j column index variable name
   * @param view_handled whether caller already handled matrix view
   * @return part of kernel with code for this and nested expressions
   */
  inline kernel_parts get_kernel_parts(
      std::set<const operation_cl_base*>& generated, name_generator& name_gen,
      const std::string& i, const std::string& j, bool view_handled) const {
    kernel_parts res{};
    if (generated.count(this) == 0) {
      this->var_name = name_gen.generate();
      generated.insert(this);

      if (PassZero && internal::matvec_mul_opt<T_no_ref>::is_possible) {
        res = internal::matvec_mul_opt<T_no_ref>::get_kernel_parts(
            this->template get_arg<0>(), generated, name_gen, i,
            var_name + "_j");
      } else {
        res = this->template get_arg<0>().get_kernel_parts(
            generated, name_gen, i, var_name + "_j", view_handled || PassZero);
      }
      kernel_parts my_part
          = generate(i, j, view_handled, this->template get_arg<0>().var_name);
      res += my_part;
      res.body = res.body_prefix + res.body;
      res.body_prefix = "";
    }
    return res;
  }

  /**
   * Generates kernel code for this expression.
   * @param i row index variable name
   * @param j column index variable name
   * @param view_handled whether whether caller already handled matrix view
   * @param var_name_arg name of the variable in kernel that holds argument to
   * this expression
   * @return part of kernel with code for this expression
   */
  inline kernel_parts generate(const std::string& i, const std::string& j,
                               const bool view_handled,
                               const std::string& var_name_arg) const {
    kernel_parts res;
    res.body_prefix
        = type_str<Scalar>() + " " + var_name + " = " + init_ + ";\n";
    if (PassZero) {
      res.body_prefix += "int " + var_name + "_start = contains_nonzero("
                         + var_name + "_view, LOWER) ? 0 : " + i + ";\n";
      if (internal::matvec_mul_opt<T_no_ref>::is_possible) {
        res.body_prefix += "int " + var_name + "_end_temp = contains_nonzero("
                           + var_name + "_view, UPPER) ? " + var_name
                           + "_cols : min(" + var_name + "_cols, " + i
                           + " + 1);\n";
        res.body_prefix += "int " + var_name + "_end = contains_nonzero("
                           + var_name + "_vec_view, UPPER) ? " + var_name
                           + "_end_temp : min(1, " + var_name + "_end_temp);\n";
      } else {
        res.body_prefix += "int " + var_name + "_end = contains_nonzero("
                           + var_name + "_view, UPPER) ? " + var_name
                           + "_cols : min(" + var_name + "_cols, " + i
                           + " + 1);\n";
      }
      res.body_prefix += "for(int " + var_name + "_j = " + var_name + "_start; "
                         + var_name + "_j < " + var_name + "_end; " + var_name
                         + "_j++){\n";
    } else {
      res.body_prefix += "for(int " + var_name + "_j = 0; " + var_name + "_j < "
                         + var_name + "_cols; " + var_name + "_j++){\n";
    }
    res.body += var_name + " = " + operation::generate(var_name, var_name_arg)
                + ";\n}\n";
    res.args = "int " + var_name + "_view, int " + var_name + "_cols, ";
    if (PassZero && internal::matvec_mul_opt<T_no_ref>::is_possible) {
      res.args += "int " + var_name + "_vec_view, ";
    }
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
  inline void set_args(std::set<const operation_cl_base*>& generated,
                       cl::Kernel& kernel, int& arg_num) const {
    if (generated.count(this) == 0) {
      generated.insert(this);
      this->template get_arg<0>().set_args(generated, kernel, arg_num);
      kernel.setArg(arg_num++, this->template get_arg<0>().view());
      kernel.setArg(arg_num++, this->template get_arg<0>().cols());
      if (PassZero && internal::matvec_mul_opt<T>::is_possible) {
        kernel.setArg(arg_num++, internal::matvec_mul_opt<T_no_ref>::view(
                                     this->template get_arg<0>()));
      }
    }
  }

  /**
   * Number of columns of a matrix that would be the result of evaluating this
   * expression.
   * @return number of columns
   */
  inline int cols() const { return 1; }

  /**
   * Determine indices of extreme sub- and superdiagonals written.
   * @return pair of indices - bottom and top diagonal
   */
  inline std::pair<int, int> extreme_diagonals() const {
    return {-rows() + 1, cols() - 1};
  }
};  // namespace math

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
  using base = rowwise_reduction<rowwise_sum_<T>, T, sum_op, true>;
  using base::arguments_;

 public:
  explicit rowwise_sum_(T&& a) : base(std::forward<T>(a), "0") {}

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& arg_copy = this->template get_arg<0>().deep_copy();
    return rowwise_sum_<std::remove_reference_t<decltype(arg_copy)>>(
        std::move(arg_copy));
  }
};

/**
 * Rowwise sum reduction of a kernel generator expression.
 * @tparam T type of input expression
 * @param a expression to reduce
 * @return sum
 */
template <typename T,
          typename = require_all_valid_expressions_and_none_scalar_t<T>>
inline auto rowwise_sum(T&& a) {
  auto&& arg_copy = as_operation_cl(std::forward<T>(a)).deep_copy();
  return rowwise_sum_<std::remove_reference_t<decltype(arg_copy)>>(
      std::move(arg_copy));
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
  using op = max_op<typename std::remove_reference_t<T>::Scalar>;
  using base = rowwise_reduction<rowwise_max_<T>, T, op, false>;
  using base::arguments_;

 public:
  explicit rowwise_max_(T&& a) : base(std::forward<T>(a), op::init()) {}
  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& arg_copy = this->template get_arg<0>().deep_copy();
    return rowwise_max_<std::remove_reference_t<decltype(arg_copy)>>(
        std::move(arg_copy));
  }
};

/**
 * Rowwise max reduction of a kernel generator expression.
 * @tparam T type of input expression
 * @param a expression to reduce
 * @return max
 */
template <typename T,
          typename = require_all_valid_expressions_and_none_scalar_t<T>>
inline auto rowwise_max(T&& a) {
  auto&& arg_copy = as_operation_cl(std::forward<T>(a)).deep_copy();
  return rowwise_max_<std::remove_reference_t<decltype(arg_copy)>>(
      std::move(arg_copy));
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
  using op = min_op<typename std::remove_reference_t<T>::Scalar>;
  using base = rowwise_reduction<rowwise_min_<T>, T, op, false>;
  using base::arguments_;

 public:
  explicit rowwise_min_(T&& a) : base(std::forward<T>(a), op::init()) {}
  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& arg_copy = this->template get_arg<0>().deep_copy();
    return rowwise_min_<std::remove_reference_t<decltype(arg_copy)>>(
        std::move(arg_copy));
  }
};

/**
 * Min reduction of a kernel generator expression.
 * @tparam T type of input expression
 * @param a expression to reduce
 * @return min
 */
template <typename T,
          typename = require_all_valid_expressions_and_none_scalar_t<T>>
inline auto rowwise_min(T&& a) {
  auto&& arg_copy = as_operation_cl(std::forward<T>(a)).deep_copy();
  return rowwise_min_<std::remove_reference_t<decltype(arg_copy)>>(
      std::move(arg_copy));
}
/** @}*/
}  // namespace math
}  // namespace stan

#endif
#endif
