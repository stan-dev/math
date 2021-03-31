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
   * @param condition condition expression
   * @param then then expression
   * @param els else expression
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

/**
 * Represents a selection operation in kernel generator expressions. This is
 * element wise ternary operator <code>condition ? then : els</code>, also
 * equivalent to Eigen's \c .select().
 * @tparam Derived derived type
 * @tparam T_condition type of condition
 * @tparam T_then type of then expression
 * @tparam T_else type of else expression
 */
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
   * @param condition condition expression
   * @param then then expression
   * @param els else expression
   */
  opencl_code_impl(names_tuple names, T_arguments&&... arguments)
      : base(std::forward<T_arguments>(arguments)...), names_(names) {
    /*for_each([this](const auto& el){
      const auto& first = this.get_arg<0>();
      if (first.rows() != base::dynamic && el.rows() != base::dynamic) {
        check_size_match("select", "Rows of ", "first argument", first.rows(),
                         "rows of ", "argument", el.rows());
      }
      if (first.cols() != base::dynamic && el.cols() != base::dynamic) {
        check_size_match("select", "Columns of ", "first argument",
    first.cols(), "columns of ", "argument", el.cols());
      }
    }, this->arguments_);*/
  }

  //  /**
  //   * Creates a deep copy of this expression.
  //   * @return copy of \c *this
  //   */
  //  inline auto deep_copy() const {
  //    return index_apply<sizeof...(T_arguments)>([this](auto... Is) {
  //      auto args_copy
  //          = std::make_tuple(this->template get_argument<Is>().deep_copy());
  //      return opencl_code_impl<
  //          Code,
  //          std::remove_reference_t<decltype(std::get<Is>(args_copy))>...>(
  //          std::move(std::get<Is>(args_copy))...);
  //    });
  //  }

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

template <const char* Code, typename... T_arguments>
class opencl_code_ : public operation_cl_base {
 public:
  std::shared_ptr<opencl_code_impl<Code, T_arguments...>> impl_;
  std::string& var_name_;
  using names_tuple
      = std::tuple<typename std::pair<const char*, T_arguments>::first_type...>;
  using Deriv = opencl_code_impl<Code, T_arguments...>;

  opencl_code_(const names_tuple& names, T_arguments&&... arguments)
      : impl_(std::make_shared<opencl_code_impl<Code, T_arguments...>>(
          names, std::forward<T_arguments>(arguments)...)),
        var_name_(impl_->var_name_) {}

  opencl_code_(const opencl_code_<Code, T_arguments...>& other)
      : impl_(other.impl_), var_name_(impl_->var_name_) {}

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
  auto set_args(std::map<const void*, const char*>& generated,
                std::map<const void*, const char*>& generated_all,
                cl::Kernel& kernel, int& arg_num) const {
    return impl_->set_args(generated, generated_all, kernel, arg_num);
  }
  auto add_read_event(cl::Event& e) const { return impl_->add_read_event(e); }
  auto get_write_events(std::vector<cl::Event>& events) const {
    return impl_->get_write_events(events);
  }
  template <std::enable_if_t<(sizeof...(T_arguments) > 0)>* = nullptr>
  auto rows() const {
    return impl_->rows();
  }
  template <std::enable_if_t<(sizeof...(T_arguments) == 0)>* = nullptr>
  auto rows() const {
    return -1;
  }
  template <std::enable_if_t<(sizeof...(T_arguments) > 0)>* = nullptr>
  auto cols() const {
    return impl_->cols();
  }
  template <std::enable_if_t<(sizeof...(T_arguments) == 0)>* = nullptr>
  auto cols() const {
    return -1;
  }
  template <std::enable_if_t<(sizeof...(T_arguments) > 0)>* = nullptr>
  auto thread_rows() const {
    return impl_->thread_rows();
  }
  template <std::enable_if_t<(sizeof...(T_arguments) == 0)>* = nullptr>
  auto thread_rows() const {
    return -1;
  }
  template <std::enable_if_t<(sizeof...(T_arguments) > 0)>* = nullptr>
  auto thread_cols() const {
    return impl_->thread_cols();
  }
  template <std::enable_if_t<(sizeof...(T_arguments) == 0)>* = nullptr>
  auto thread_cols() const {
    return -1;
  }
  auto extreme_diagonals() const { return impl_->extreme_diagonals(); }
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
      auto args_copy = std::make_tuple(
          this->impl->template get_argument<Is>().deep_copy()...);
      return opencl_code_<
          Code, std::remove_reference_t<decltype(std::get<Is>(args_copy))>...>(
          this->impl_->names, std::move(std::get<Is>(args_copy))...);
    });
  }

  template <typename T_scalar>
  inline auto output(const char* var_name) const {
    return opencl_code_output<opencl_code_<Code, T_arguments...>, T_scalar>(
        *this, var_name);
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
template <const char* Code, typename... T_arguments,
          require_all_kernel_expressions_t<T_arguments...>* = nullptr>
inline opencl_code_<Code, as_operation_cl_t<T_arguments>...> opencl_code(
    std::tuple<typename std::pair<const char*, T_arguments>::first_type...>
        names,
    T_arguments&&... arguments) {
  return {names, as_operation_cl(std::forward<T_arguments>(arguments))...};
}

/** @}*/
}  // namespace math
}  // namespace stan
#endif
#endif
