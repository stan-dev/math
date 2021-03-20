#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_ROWWISE_REDUCTION_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_ROWWISE_REDUCTION_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/broadcast.hpp>
#include <stan/math/opencl/kernel_generator/binary_operation.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <map>
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
      const Arg& a, std::map<const void*, const char*>& generated,
      std::map<const void*, const char*>& generated_all,
      name_generator& name_gen, const std::string& row_index_name,
      const std::string& col_index_name) {
    return {};
  }
};

template <typename Mat, typename VecT>
struct matvec_mul_opt<elt_multiply_<Mat, broadcast_<VecT, true, false>>> {
  // if the argument of rowwise reduction is multiplication with a broadcast
  // vector we can do the optimization
  enum { is_possible = 1 };
  using Arg = elt_multiply_<Mat, broadcast_<VecT, true, false>>;

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
   * @param[in,out] generated map from (pointer to) already generated local
   * operations to variable names
   * @param[in,out] generated_all map from (pointer to) already generated all
   * operations to variable names
   * @param name_gen name generator for this kernel
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @return part of kernel with code for this and nested expressions
   */
  static kernel_parts get_kernel_parts(
      const Arg& mul, std::map<const void*, const char*>& generated,
      std::map<const void*, const char*>& generated_all,
      name_generator& name_gen, const std::string& row_index_name,
      const std::string& col_index_name) {
    kernel_parts res{};
    if (generated.count(&mul) == 0) {
      mul.var_name_ = name_gen.generate();
      generated[&mul] = "";

      const auto& matrix = mul.template get_arg<0>();
      const auto& broadcast = mul.template get_arg<1>();
      res = matrix.get_kernel_parts(generated, generated_all, name_gen,
                                    row_index_name, col_index_name, true);
      if (generated.count(&broadcast) == 0) {
        broadcast.var_name_ = name_gen.generate();
        generated[&broadcast] = "";

        const auto& vec = broadcast.template get_arg<0>();
        std::string row_index_name_bc = row_index_name;
        std::string col_index_name_bc = col_index_name;
        broadcast.modify_argument_indices(row_index_name_bc, col_index_name_bc);
        res += vec.get_kernel_parts(generated, generated_all, name_gen,
                                    row_index_name_bc, col_index_name_bc, true);
        res += broadcast.generate(row_index_name, col_index_name, true,
                                  vec.var_name_);
      }
      res += mul.generate(row_index_name, col_index_name, true,
                          matrix.var_name_, broadcast.var_name_);
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
  using base::var_name_;

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
   * @param[in,out] generated map from (pointer to) already generated local
   * operations to variable names
   * @param[in,out] generated_all map from (pointer to) already generated all
   * operations to variable names
   * @param name_gen name generator for this kernel
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @param view_handled whether caller already handled matrix view
   * @return part of kernel with code for this and nested expressions
   */
  inline kernel_parts get_kernel_parts(
      std::map<const void*, const char*>& generated,
      std::map<const void*, const char*>& generated_all,
      name_generator& name_gen, const std::string& row_index_name,
      const std::string& col_index_name, bool view_handled) const {
    kernel_parts res{};
    if (generated.count(this) == 0) {
      this->var_name_ = name_gen.generate();
      generated[this] = "";

      std::map<const void*, const char*> generated2;
      if (PassZero && internal::matvec_mul_opt<T_no_ref>::is_possible) {
        res = internal::matvec_mul_opt<T_no_ref>::get_kernel_parts(
            this->template get_arg<0>(), generated2, generated_all, name_gen,
            row_index_name, var_name_ + "_j");
      } else {
        res = this->template get_arg<0>().get_kernel_parts(
            generated2, generated_all, name_gen, row_index_name,
            var_name_ + "_j", view_handled || PassZero);
      }
      kernel_parts my_part
          = generate(row_index_name, col_index_name, view_handled,
                     this->template get_arg<0>().var_name_);
      res += my_part;
      res.body = res.body_prefix + res.body;
      res.body_prefix = "";
    }
    return res;
  }

  /**
   * Generates kernel code for this expression.
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @param view_handled whether whether caller already handled matrix view
   * @param var_name_arg name of the variable in kernel that holds argument to
   * this expression
   * @return part of kernel with code for this expression
   */
  inline kernel_parts generate(const std::string& row_index_name,
                               const std::string& col_index_name,
                               const bool view_handled,
                               const std::string& var_name_arg) const {
    kernel_parts res;
    res.body_prefix
        = type_str<Scalar>() + " " + var_name_ + " = " + init_ + ";\n";
    if (PassZero) {
      res.body_prefix += "int " + var_name_ + "_start = contains_nonzero("
                         + var_name_ + "_view, LOWER) ? 0 : " + row_index_name
                         + ";\n";
      if (internal::matvec_mul_opt<T_no_ref>::is_possible) {
        res.body_prefix += "int " + var_name_ + "_end_temp = contains_nonzero("
                           + var_name_ + "_view, UPPER) ? " + var_name_
                           + "_cols : min(" + var_name_ + "_cols, "
                           + row_index_name + " + 1);\n";
        res.body_prefix += "int " + var_name_ + "_end = contains_nonzero("
                           + var_name_ + "_vec_view, UPPER) ? " + var_name_
                           + "_end_temp : min(1, " + var_name_
                           + "_end_temp);\n";
      } else {
        res.body_prefix += "int " + var_name_ + "_end = contains_nonzero("
                           + var_name_ + "_view, UPPER) ? " + var_name_
                           + "_cols : min(" + var_name_ + "_cols, "
                           + row_index_name + " + 1);\n";
      }
      res.body_prefix += "for(int " + var_name_ + "_j = " + var_name_
                         + "_start; " + var_name_ + "_j < " + var_name_
                         + "_end; " + var_name_ + "_j++){\n";
    } else {
      res.body_prefix += "for(int " + var_name_ + "_j = 0; " + var_name_
                         + "_j < " + var_name_ + "_cols; " + var_name_
                         + "_j++){\n";
    }
    res.body += var_name_ + " = " + operation::generate(var_name_, var_name_arg)
                + ";\n}\n";
    res.args = "int " + var_name_ + "_view, int " + var_name_ + "_cols, ";
    if (PassZero && internal::matvec_mul_opt<T_no_ref>::is_possible) {
      res.args += "int " + var_name_ + "_vec_view, ";
    }
    return res;
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
 * @param a the expression to reduce
 * @return sum
 */
template <typename T,
          typename = require_all_kernel_expressions_and_none_scalar_t<T>>
inline auto rowwise_sum(T&& a) {
  auto&& arg_copy = as_operation_cl(std::forward<T>(a)).deep_copy();
  return rowwise_sum_<std::remove_reference_t<decltype(arg_copy)>>(
      std::move(arg_copy));
}

/**
 * Operation for product reduction.
 */
struct prod_op {
  /**
   * Generates prod reduction kernel code.
   * @param a first variable
   * @param b second variable
   * @return reduction code
   */
  inline static std::string generate(const std::string& a,
                                     const std::string& b) {
    return a + " * " + b;
  }
};

/**
 * Represents rowwise product reduction in kernel generator expressions.
 * @tparam T type of expression
 */
template <typename T>
class rowwise_prod_
    : public rowwise_reduction<rowwise_prod_<T>, T, prod_op, false> {
  using base = rowwise_reduction<rowwise_prod_<T>, T, prod_op, false>;
  using base::arguments_;

 public:
  explicit rowwise_prod_(T&& a) : base(std::forward<T>(a), "1") {}

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& arg_copy = this->template get_arg<0>().deep_copy();
    return rowwise_prod_<std::remove_reference_t<decltype(arg_copy)>>(
        std::move(arg_copy));
  }
};

/**
 * Rowwise product reduction of a kernel generator expression.
 * @tparam T type of input expression
 * @param a the expression to reduce
 * @return prod
 */
template <typename T,
          typename = require_all_kernel_expressions_and_none_scalar_t<T>>
inline auto rowwise_prod(T&& a) {
  auto&& arg_copy = as_operation_cl(std::forward<T>(a)).deep_copy();
  return rowwise_prod_<std::remove_reference_t<decltype(arg_copy)>>(
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
 * @param a the expression to reduce
 * @return max
 */
template <typename T,
          typename = require_all_kernel_expressions_and_none_scalar_t<T>>
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
 * @param a the expression to reduce
 * @return min
 */
template <typename T,
          typename = require_all_kernel_expressions_and_none_scalar_t<T>>
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
