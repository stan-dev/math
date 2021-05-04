#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_OPERATION_CL_LHS_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_OPERATION_CL_LHS_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <string>
#include <map>
#include <array>
#include <numeric>
#include <vector>

namespace stan {
namespace math {

/** \addtogroup opencl_kernel_generator
 *  @{
 */
/**
 * Base for all kernel generator operations that can be used on left hand side
 * of an expression.
 * @tparam Derived derived type
 * @tparam Scalar scalar type of the result
 * @tparam Args types of arguments to this operation
 */
template <typename Derived, typename Scalar, typename... Args>
class operation_cl_lhs : public operation_cl<Derived, Scalar, Args...>,
                         public operation_cl_lhs_base {
 protected:
  using base = operation_cl<Derived, Scalar, Args...>;
  static constexpr int N = sizeof...(Args);

 public:
  using base::derived;
  using base::operation_cl;

  /**
   * Generates kernel code for this expression if it appears on the left hand
   * side of an assignment.
   * @param[in,out] generated map from (pointer to) already generated local
   * operations to variable names
   * @param[in,out] generated_all map from (pointer to) already generated all
   * operations to variable names
   * @param name_gen name generator for this kernel
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @return part of kernel with code for this expressions
   */
  inline kernel_parts get_kernel_parts_lhs(
      std::map<const void*, const char*>& generated,
      std::map<const void*, const char*>& generated_all,
      name_generator& name_gen, const std::string& row_index_name,
      const std::string& col_index_name) const {
    if (generated.count(this) == 0) {
      generated[this] = "";
      this->var_name_ = name_gen.generate();
    }
    std::string row_index_name_arg = row_index_name;
    std::string col_index_name_arg = col_index_name;
    derived().modify_argument_indices(row_index_name_arg, col_index_name_arg);
    std::array<kernel_parts, N> args_parts = index_apply<N>([&](auto... Is) {
      std::map<const void*, const char*> generated2;
      return std::array<kernel_parts, N>{
          this->template get_arg<Is>().get_kernel_parts_lhs(
              &Derived::modify_argument_indices
                      == &operation_cl<Derived, Scalar,
                                       Args...>::modify_argument_indices
                  ? generated
                  : generated2,
              generated_all, name_gen, row_index_name_arg,
              col_index_name_arg)...};
    });
    kernel_parts res
        = std::accumulate(args_parts.begin(), args_parts.end(), kernel_parts{});
    kernel_parts my_part = index_apply<N>([&](auto... Is) {
      return this->derived().generate_lhs(
          row_index_name, col_index_name,
          this->template get_arg<Is>().var_name_...);
    });
    res += my_part;
    if (generated_all.count(this) == 0) {
      generated_all[this] = "";
    } else {
      res.args = "";
    }
    return res;
  }

  /**
   * Generates kernel code for this and nested expressions if this expression
   * appears on the left hand side of an assignment.
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @param var_name_arg name of the variable in kernel that holds argument to
   * this expression
   * @return part of kernel with code for this expression
   */
  inline kernel_parts generate_lhs(const std::string& row_index_name,
                                   const std::string& col_index_name,
                                   const std::string& var_name_arg) const {
    return {};
  }

  /**
   * Evaluates an expression and assigns it to this.
   * @tparam T_expression type of expression
   * @param rhs input expression
   */
  template <typename T_expression,
            typename
            = require_all_kernel_expressions_and_none_scalar_t<T_expression>>
  Derived& operator=(T_expression&& rhs) {
    auto expression
        = as_operation_cl(std::forward<T_expression>(rhs)).derived();
    int this_rows = derived().rows();
    int this_cols = derived().cols();
    if (this_rows == expression.rows() && this_cols == expression.cols()
        && this_rows * this_cols == 0) {
      return derived();
    }
    expression.evaluate_into(derived());
    return derived();
  }
  // Copy assignment delegates to general assignment operator. If we didn't
  // implement this, we would get ambiguities in overload resolution with
  // implicitly generated one
  inline const operation_cl_lhs<Derived, Scalar, Args...>& operator=(
      const operation_cl_lhs<Derived, Scalar, Args...>& rhs) const {
    return operator=<const operation_cl_lhs<Derived, Scalar, Args...>&>(rhs);
  }

  /**
   * Sets the view of the underlying matrix depending on which of its parts are
   * written to.
   * @param bottom_diagonal Index of the top sub- or super- diagonal written
   * with nonzero elements.
   * @param top_diagonal Index of the top sub- or super- diagonal written with
   * nonzero elements.
   * @param bottom_zero_diagonal Index of the top sub- or super- diagonal
   * written with zeros if it ie more extreme than \c bottom_diagonal. Otherwise
   * it should be set to equal value as \c bottom_diagonal.
   * @param top_zero_diagonal Index of the top sub- or super- diagonal written
   * with zeros if it ie more extreme than \c top_diagonal. Otherwise it should
   * be set to equal value as \c top_diagonal.
   */
  inline void set_view(int bottom_diagonal, int top_diagonal,
                       int bottom_zero_diagonal, int top_zero_diagonal) const {
    index_apply<N>([&](auto... Is) {
      static_cast<void>(std::initializer_list<int>{
          (this->template get_arg<Is>().set_view(bottom_diagonal, top_diagonal,
                                                 bottom_zero_diagonal,
                                                 top_zero_diagonal),
           0)...});
    });
  }

  /**
   * Sets the dimensions of the underlying expressions if possible. If not
   * checks whether they have correct dimensions.
   * @param rows desired number of rows
   * @param cols desired number of columns
   * @throws std::invalid_argument desired dimensions do not match with
   * dimensions of underlying expression that can not be resized.
   */
  inline void check_assign_dimensions(int rows, int cols) const {
    index_apply<N>([&](auto... Is) {
      static_cast<void>(std::initializer_list<int>{
          (this->template get_arg<Is>().check_assign_dimensions(rows, cols),
           0)...});
    });
  }

  /**
   * Adds write event to any matrices used by nested expressions.
   * @param e the event to add
   */
  inline void add_write_event(cl::Event& e) const {
    index_apply<N>([&](auto... Is) {
      static_cast<void>(std::initializer_list<int>{
          (this->template get_arg<Is>().add_write_event(e), 0)...});
    });
  }

  /**
   * Adds all read and write events on any matrices used by nested expressions
   * to a list and clears them from those matrices.
   * @param[out] events List of all events.
   */
  inline void get_clear_read_write_events(
      std::vector<cl::Event>& events) const {
    index_apply<N>([&](auto... Is) {
      static_cast<void>(std::initializer_list<int>{
          (this->template get_arg<Is>().get_clear_read_write_events(events),
           0)...});
    });
  }
};
/** @}*/
}  // namespace math
}  // namespace stan

#endif
#endif
