#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_OPENCL_CODE_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_OPENCL_CODE_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/common_return_scalar.hpp>
#include <stan/math/prim/functor/for_each.hpp>
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
 * Represents output variable of custom code in kernel generator expressions.
 * @tparam T_code intantiation of `opencl_code_` class template this is output
 * of
 * @tparam T_scalar scalar type of the output variable
 */
template <typename T_code, typename T_scalar>
class opencl_code_output
    : public operation_cl<opencl_code_output<T_code, T_scalar>, T_scalar,
                          T_code> {
 public:
  using Scalar = T_scalar;
  using base
      = operation_cl<opencl_code_output<T_code, T_scalar>, T_scalar, T_code>;
  using base::var_name_;
  const char* custom_var_name_;

  /**
   * Constructor
   * @param code code object
   * @param custom_var_name variable name from the custom code
   */
  opencl_code_output(T_code code, const char* custom_var_name)
      : base(std::move(code)), custom_var_name_(custom_var_name) {}

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& code_copy = this->template get_arg<0>().deep_copy();
    return opencl_code_output<std::remove_reference_t<decltype(code_copy)>,
                              T_scalar>(std::move(code_copy), custom_var_name_);
  }

  /**
   * Generates kernel code for this operation.
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @param view_handled whether whether caller already handled matrix view
   * @param dummy_var_name_code variable name of the code operation (unused)
   * @return part of kernel with code for this expression
   */
  inline kernel_parts generate(const std::string& row_index_name,
                               const std::string& col_index_name,
                               const bool view_handled,
                               const std::string& dummy_var_name_code) const {
    kernel_parts res{};
    res.body = type_str<Scalar>() + " " + var_name_ + " = " + custom_var_name_
               + ";\n";
    return res;
  }

  /**
   * Determine indices of extreme sub- and superdiagonals written.
   * @return pair of indices - bottom and top diagonal
   */
  inline std::pair<int, int> extreme_diagonals() const {
    return {-this->rows() + 1, this->cols() - 1};
  }
};

namespace internal {
template <const char* Code, typename... T_arguments>
class opencl_code_impl
    : public operation_cl<opencl_code_impl<Code, T_arguments...>, double,
                          T_arguments...> {
 public:
  //  using Scalar = double;
  using base = operation_cl<opencl_code_impl<Code, T_arguments...>, double,
                            T_arguments...>;
  using names_tuple
      = std::tuple<typename std::pair<const char*, T_arguments>::first_type...>;
  using base::var_name_;
  names_tuple names_;

  /**
   * Constructor
   * @param names tuple of names of the input variables (that hold values of
   * argument expressions)
   * @param arguments arguments to this operation
   */
  explicit opencl_code_impl(names_tuple names, T_arguments&&... arguments)
      : base(std::forward<T_arguments>(arguments)...), names_(names) {}

  /**
   * Generates kernel code for this (select) operation.
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @param view_handled whether whether caller already handled matrix
   * view
   * @param var_names variable names of the argument expressions
   * @return part of kernel with code for this expression
   */
  inline kernel_parts generate(
      const std::string& row_index_name, const std::string& col_index_name,
      const bool view_handled,
      std::tuple_element_t<
          0, std::pair<const std::string&, T_arguments>>... var_names) const {
    return index_apply<sizeof...(T_arguments)>([this](auto... Is) {
      kernel_parts res{};
      std::array<std::string, sizeof...(T_arguments)> input_renames{
          (type_str<scalar_type_t<decltype(this->template get_arg<Is>())>>()
           + " " + std::get<Is>(names_) + " = "
           + this->template get_arg<Is>().var_name_ + ";\n")...};
      res.body = std::accumulate(input_renames.begin(), input_renames.end(),
                                 std::string())
                 + Code;
      return res;
    });
  }
};
}  // namespace internal

/**
 * Represents custom code in kernel generator expressions.
 * @tparam Code custom code
 * @tparam T_arguments types of argument expressions
 */
template <const char* Code, typename... T_arguments>
class opencl_code_ : public operation_cl_base {
 public:
  std::shared_ptr<internal::opencl_code_impl<Code, T_arguments...>> impl_;
  std::string& var_name_;
  using names_tuple
      = std::tuple<typename std::pair<const char*, T_arguments>::first_type...>;
  using Deriv = internal::opencl_code_impl<Code, T_arguments...>;

  /**
   * Constructor
   * @param names tuple of names of the input variables (that hold values of
   * argument expressions)
   * @param arguments arguments to this operation
   */
  explicit opencl_code_(const names_tuple& names, T_arguments&&... arguments)
      : impl_(
            std::make_shared<internal::opencl_code_impl<Code, T_arguments...>>(
                names, std::forward<T_arguments>(arguments)...)),
        var_name_(impl_->var_name_) {}

  /**
   * Copy constructor.
   * @param other object to copy
   */
  opencl_code_(const opencl_code_<Code, T_arguments...>& other)
      : impl_(other.impl_), var_name_(impl_->var_name_) {}

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
  auto get_kernel_parts(std::map<const void*, const char*>& generated,
                        std::map<const void*, const char*>& generated_all,
                        name_generator& name_gen,
                        const std::string& row_index_name,
                        const std::string& col_index_name,
                        bool view_handled) const {
    return impl_->get_kernel_parts(generated, generated_all, name_gen,
                                   row_index_name, col_index_name,
                                   view_handled);
  }

  /**
   * Sets kernel arguments for nested expressions.
   * @param[in,out] generated map from (pointer to) already generated local
   * operations to variable names
   * @param[in,out] generated_all map from (pointer to) already generated all
   * operations to variable names
   * @param kernel kernel to set arguments on
   * @param[in,out] arg_num consecutive number of the first argument to set.
   * This is incremented for each argument set by this function.
   */
  auto set_args(std::map<const void*, const char*>& generated,
                std::map<const void*, const char*>& generated_all,
                cl::Kernel& kernel, int& arg_num) const {
    return impl_->set_args(generated, generated_all, kernel, arg_num);
  }

  /**
   * Adds read event to any matrices used by nested expressions.
   * @param e the event to add
   */
  auto add_read_event(cl::Event& e) const { return impl_->add_read_event(e); }

  /**
   * Adds all write events on any matrices used by nested expressions to a list.
   * @param[out] events List of all events.
   */
  auto get_write_events(std::vector<cl::Event>& events) const {
    return impl_->get_write_events(events);
  }

  /**
   * Number of rows of a matrix that would be the result of evaluating this
   * expression. Some subclasses may need to override this.
   * @return number of rows
   */
  template <int N = sizeof...(T_arguments),
            std::enable_if_t<(N > 0)>* = nullptr>
  auto rows() const {
    return impl_->rows();
  }
  template <int N = sizeof...(T_arguments),
            std::enable_if_t<(N == 0)>* = nullptr>
  auto rows() const {
    return -1;
  }

  /**
   * Number of columns of a matrix that would be the result of evaluating this
   * expression. Some subclasses may need to override this.
   * @return number of columns
   */
  template <int N = sizeof...(T_arguments),
            std::enable_if_t<(N > 0)>* = nullptr>
  auto cols() const {
    return impl_->cols();
  }
  template <int N = sizeof...(T_arguments),
            std::enable_if_t<(N == 0)>* = nullptr>
  auto cols() const {
    return -1;
  }

  /**
   * Number of rows threads need to be launched for. For most expressions this
   * equals number of rows of the result.
   * @return number of rows
   */
  template <int N = sizeof...(T_arguments),
            std::enable_if_t<(N > 0)>* = nullptr>
  auto thread_rows() const {
    return impl_->thread_rows();
  }
  template <int N = sizeof...(T_arguments),
            std::enable_if_t<(N == 0)>* = nullptr>
  auto thread_rows() const {
    return -1;
  }

  /**
   * Number of columns threads need to be launched for. For most expressions
   * this equals number of cols of the result.
   * @return number of columns
   */
  template <int N = sizeof...(T_arguments),
            std::enable_if_t<(N > 0)>* = nullptr>
  auto thread_cols() const {
    return impl_->thread_cols();
  }
  template <int N = sizeof...(T_arguments),
            std::enable_if_t<(N == 0)>* = nullptr>
  auto thread_cols() const {
    return -1;
  }

  /**
   * Determine indices of extreme sub- and superdiagonals written.
   * @return pair of indices - bottom and top diagonal
   */
  auto extreme_diagonals() const { return impl_->extreme_diagonals(); }

  /**
   * Collects data that is needed beside types to uniqly identify a kernel
   * generator expression.
   * @param[out] uids ids of unique matrix accesses
   * @param[in,out] id_map map from memory addresses to unique ids
   * @param[in,out] next_id neqt unique id to use
   */
  auto get_unique_matrix_accesses(std::vector<int>& uids,
                                  std::map<const void*, int>& id_map,
                                  int& next_id) const {
    return impl_->get_unique_matrix_accesses(uids, id_map, next_id);
  }

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    return index_apply<sizeof...(T_arguments)>([this](auto... Is) {
      auto args_copy
          = std::make_tuple(this->impl_->template get_arg<Is>().deep_copy()...);
      return opencl_code_<
          Code, std::remove_reference_t<decltype(std::get<Is>(args_copy))>...>(
          this->impl_->names_, std::move(std::get<Is>(args_copy))...);
    });
  }

  /**
   * Get object representing output variable of ccustom code.
   * @param var_name name of the variable output object represents
   */
  template <typename T_scalar>
  inline auto output(const char* var_name) const {
    return opencl_code_output<opencl_code_<Code, T_arguments...>, T_scalar>(
        *this, var_name);
  }
};

/**
 * Custom code in kernel generator expressions.
 * @tparam Code custom code
 * @tparam T_arguments types of argument expressions
 * @param names tuple of names of the input variables (that hold values of
 * argument expressions)
 * @param arguments arguments to this operation
 */
template <const char* Code, typename... T_arguments,
          require_all_kernel_expressions_t<T_arguments...>* = nullptr>
inline auto opencl_code(
    std::tuple<typename std::pair<const char*, T_arguments>::first_type...>
        names,
    T_arguments&&... arguments) {
  return opencl_code_<Code, as_operation_cl_t<T_arguments>...>(
      names, as_operation_cl(std::forward<T_arguments>(arguments))...);
}

/** @}*/
}  // namespace math
}  // namespace stan
#endif
#endif
