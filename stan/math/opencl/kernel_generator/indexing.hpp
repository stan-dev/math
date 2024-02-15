#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_indexing_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_indexing_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <map>
#include <string>
#include <utility>

namespace stan {
namespace math {

/** \addtogroup opencl_kernel_generator
 *  @{
 */

/**
 * Represents indexing of a matrix with two matrices of indices. `indexing(mat,
 * col_index, row_index)[i,j] = mat[col_index[i,j], row_index[i,j]]`
 * @tparam T_mat expression type of the matrix to index
 * @tparam T_row_index expression type of the row index
 * @tparam T_col_index expression type of the column index
 */
template <typename T_mat, typename T_row_index, typename T_col_index>
class indexing_
    : public operation_cl_lhs<indexing_<T_mat, T_row_index, T_col_index>,
                              typename std::remove_reference_t<T_mat>::Scalar,
                              T_mat, T_row_index, T_col_index> {
  static_assert(std::is_integral<value_type_t<T_row_index>>::value,
                "indexing: Row index scalar type must be an integer!");
  static_assert(std::is_integral<value_type_t<T_col_index>>::value,
                "indexing: Column index scalar type must be an integer!");

 public:
  using Scalar = typename std::remove_reference_t<T_mat>::Scalar;
  using base = operation_cl_lhs<indexing_<T_mat, T_row_index, T_col_index>,
                                Scalar, T_mat, T_row_index, T_col_index>;
  using base::var_name_;
  using base::operator=;

  /**
   * Constructor
   * @param mat expression to index
   * @param row_index row index expression
   * @param col_index column index expression
   */
  indexing_(T_mat&& mat, T_row_index&& row_index, T_col_index&& col_index)
      : base(std::forward<T_mat>(mat), std::forward<T_row_index>(row_index),
             std::forward<T_col_index>(col_index)) {
    const char* function = "indexing";
    if (col_index.rows() != base::dynamic
        && row_index.rows() != base::dynamic) {
      check_size_match(function, "Rows of ", "col_index", col_index.rows(),
                       "rows of ", "row_index", row_index.rows());
    }
    if (col_index.cols() != base::dynamic
        && row_index.cols() != base::dynamic) {
      check_size_match(function, "Columns of ", "col_index", col_index.cols(),
                       "columns of ", "row_index", row_index.cols());
    }
  }

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& mat_copy = this->template get_arg<0>().deep_copy();
    auto&& row_index_copy = this->template get_arg<1>().deep_copy();
    auto&& col_index_copy = this->template get_arg<2>().deep_copy();
    return indexing_<std::remove_reference_t<decltype(mat_copy)>,
                     std::remove_reference_t<decltype(row_index_copy)>,
                     std::remove_reference_t<decltype(col_index_copy)>>{
        std::move(mat_copy), std::move(row_index_copy),
        std::move(col_index_copy)};
  }

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
      std::unordered_map<const void*, const char*>& generated,
      std::unordered_map<const void*, const char*>& generated_all,
      name_generator& name_gen, const std::string& row_index_name,
      const std::string& col_index_name, bool view_handled) const {
    kernel_parts res{};
    if (generated.count(this) == 0) {
      generated[this] = "";

      const auto& mat = this->template get_arg<0>();
      const auto& row_index = this->template get_arg<1>();
      const auto& col_index = this->template get_arg<2>();

      kernel_parts parts_row_idx = row_index.get_kernel_parts(
          generated, generated_all, name_gen, row_index_name, col_index_name,
          view_handled);
      kernel_parts parts_col_idx = col_index.get_kernel_parts(
          generated, generated_all, name_gen, row_index_name, col_index_name,
          view_handled);
      std::unordered_map<const void*, const char*> generated2;
      kernel_parts parts_mat = mat.get_kernel_parts(
          generated2, generated_all, name_gen, row_index.var_name_,
          col_index.var_name_, false);

      res = parts_row_idx + parts_col_idx + parts_mat;
      var_name_ = mat.var_name_;
    }
    return res;
  }

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
      std::unordered_map<const void*, const char*>& generated,
      std::unordered_map<const void*, const char*>& generated_all,
      name_generator& name_gen, const std::string& row_index_name,
      const std::string& col_index_name) const {
    if (generated.count(this) == 0) {
      generated[this] = "";
    }
    const auto& mat = this->template get_arg<0>();
    const auto& row_index = this->template get_arg<1>();
    const auto& col_index = this->template get_arg<2>();

    kernel_parts parts_row_idx
        = row_index.get_kernel_parts(generated, generated_all, name_gen,
                                     row_index_name, col_index_name, false);
    kernel_parts parts_col_idx
        = col_index.get_kernel_parts(generated, generated_all, name_gen,
                                     row_index_name, col_index_name, false);
    std::unordered_map<const void*, const char*> generated2;
    kernel_parts parts_mat
        = mat.get_kernel_parts_lhs(generated2, generated_all, name_gen,
                                   row_index.var_name_, col_index.var_name_);

    kernel_parts res = parts_row_idx + parts_col_idx + parts_mat;
    var_name_ = mat.var_name_;
    return res;
  }

  /**
   * Sets kernel arguments for this expression.
   * @param[in,out] generated map from (pointer to) already generated local
   * operations to variable names
   * @param[in,out] generated_all map from (pointer to) already generated all
   * operations to variable names
   * @param kernel kernel to set arguments on
   * @param[in,out] arg_num consecutive number of the first argument to set.
   * This is incremented for each argument set by this function.
   */
  inline void set_args(
      std::unordered_map<const void*, const char*>& generated,
      std::unordered_map<const void*, const char*>& generated_all,
      cl::Kernel& kernel, int& arg_num) const {
    if (generated.count(this) == 0) {
      generated[this] = "";
      this->template get_arg<1>().set_args(generated, generated_all, kernel,
                                           arg_num);
      this->template get_arg<2>().set_args(generated, generated_all, kernel,
                                           arg_num);
      std::unordered_map<const void*, const char*> generated2;
      this->template get_arg<0>().set_args(generated2, generated_all, kernel,
                                           arg_num);
    }
  }

  /**
   * Number of rows of a matrix that would be the result of evaluating this
   * expression.
   * @return number of rows
   */
  inline int rows() const {
    return std::max(this->template get_arg<1>().rows(),
                    this->template get_arg<2>().rows());
  }

  /**
   * Number of columns of a matrix that would be the result of evaluating this
   * expression.
   * @return number of columns
   */
  inline int cols() const {
    return std::max(this->template get_arg<1>().cols(),
                    this->template get_arg<2>().cols());
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
    this->template get_arg<0>().set_view(
        std::numeric_limits<int>::min(), std::numeric_limits<int>::max(),
        std::numeric_limits<int>::min(), std::numeric_limits<int>::max());
  }

  /**
   * Determine indices of extreme sub- and superdiagonals written.
   * @return pair of indices - bottom and top diagonal
   */
  inline std::pair<int, int> extreme_diagonals() const {
    return {std::numeric_limits<int>::min(), std::numeric_limits<int>::max()};
  }

  /**
   * Checks if desired dimensions match dimensions of the indexing.
   * @param rows desired number of rows
   * @param cols desired number of columns
   * @throws std::invalid_argument desired dimensions do not match dimensions
   * of the indexing.
   */
  inline void check_assign_dimensions(int rows, int cols) const {
    check_size_match("indexing_.check_assign_dimensions", "Rows of ",
                     "indexing", this->rows(), "rows of ", "expression", rows);
    check_size_match("indexing_.check_assign_dimensions", "Columns of ",
                     "indexing", this->cols(), "columns of ", "expression",
                     cols);
  }

  /**
   * Adds write event to indexed matrix and read event to indices.
   * @param e the event to add
   */
  inline void add_write_event(cl::Event& e) const {
    this->template get_arg<1>().add_read_event(e);
    this->template get_arg<2>().add_read_event(e);
    this->template get_arg<0>().add_write_event(e);
  }

  /**
   * Adds all read and write events on any matrices used by nested expressions
   * to a list and clears them from those matrices.
   * @param[out] events List of all events.
   */
  inline void get_clear_read_write_events(
      std::vector<cl::Event>& events) const {
    this->template get_arg<0>().get_clear_read_write_events(events);
    this->template get_arg<1>().get_write_events(events);
    this->template get_arg<2>().get_write_events(events);
  }
};

/**
 * Index a kernel generator expression using two expressions for indices. The
 * result is a matrix of the same size as index matrices and with same scalar
 * type as indexed expression. `indexing(mat, col_index, row_index)[i,j] =
 * mat[col_index[i,j], row_index[i,j]]`
 *
 * If a matrix is both indexed and the result of the kernel (such as in
 * `indexing(a, b, c) = indexing(a, d, e);`), the result can be wrong due to
 * aliasing. In this case the expression should be evaluated in a temporary by
 * doing `indexing(a, b, c) = indexing(a, d, e).eval();`. This is not necessary
 * if both indexings use the same indices or index no common elements of the
 * matrix.
 *
 * If indexing is assigned to and some element is indexed multiple times it can
 * end with either of the assigned values due to a data race.
 *
 * @tparam T_mat type of indexed expression
 * @tparam T_row_index type of row index
 * @tparam T_col_index type of column index
 * @param mat indexed matrix
 * @param row_index row index
 * @param col_index column index
 * @return indexing expression
 */
template <typename T_mat, typename T_row_index, typename T_col_index,
          require_all_kernel_expressions_t<T_mat, T_row_index,
                                           T_col_index>* = nullptr>
inline auto indexing(T_mat&& mat, T_row_index&& row_index,
                     T_col_index&& col_index) {
  auto&& mat_operation = as_operation_cl(std::forward<T_mat>(mat)).deep_copy();
  return indexing_<std::remove_reference_t<decltype(mat_operation)>,
                   as_operation_cl_t<T_row_index>,
                   as_operation_cl_t<T_col_index>>(
      std::move(mat_operation),
      as_operation_cl(std::forward<T_row_index>(row_index)),
      as_operation_cl(std::forward<T_col_index>(col_index)));
}

/** @}*/
}  // namespace math
}  // namespace stan

#endif
#endif
