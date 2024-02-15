#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_LOAD_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_LOAD_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernel_generator/assignment_ops.hpp>

#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl_lhs.hpp>
#include <type_traits>
#include <string>
#include <utility>
#include <map>
#include <vector>

namespace stan {
namespace math {

/** \addtogroup opencl_kernel_generator
 *  @{
 */
/**
 * Represents an access to a \c matrix_cl in kernel generator expressions
 * @tparam T \c matrix_cl
 * @tparam AssignOp tells higher level operations whether the final operation
 * should be an assignment or a type of compound assignment.
 */
template <typename T, assign_op_cl AssignOp = assign_op_cl::equals>
class load_
    : public operation_cl_lhs<load_<T, AssignOp>,
                              typename std::remove_reference_t<T>::type> {
 protected:
  T a_;

 public:
  static constexpr assign_op_cl assignment_op = AssignOp;
  using Scalar = typename std::remove_reference_t<T>::type;
  using base = operation_cl<load_<T, AssignOp>, Scalar>;
  using base::var_name_;
  static_assert(disjunction<is_matrix_cl<T>, is_arena_matrix_cl<T>>::value,
                "load_: argument a must be a matrix_cl<T>!");
  static_assert(
      std::is_arithmetic<Scalar>::value,
      "load_: T in \"matrix_cl<T> a\" argument must be an arithmetic type!");

  /**
   * Constructor
   * @param a \c matrix_cl
   */
  explicit load_(T&& a) : a_(std::forward<T>(a)) {}

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline load_<T&, AssignOp> deep_copy() & { return load_<T&, AssignOp>(a_); }
  inline load_<const T&, AssignOp> deep_copy() const& {
    return load_<const T&, AssignOp>(a_);
  }
  inline load_<T, AssignOp> deep_copy() && {
    return load_<T, AssignOp>(std::forward<T>(a_));
  }

  /**
   * Generates kernel code for this expression.
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
    if (generated.count(&a_) == 0) {
      if (generated_all.count(&a_) == 0) {
        this->var_name_ = name_gen.generate();
        generated_all[&a_] = generated[&a_] = this->var_name_.c_str();
        res.body = generate_body(row_index_name, col_index_name, view_handled,
                                 this->var_name_.c_str());
        res.args = "__global " + type_str<Scalar>() + "* " + var_name_
                   + "_global, int " + var_name_ + "_rows, int " + var_name_
                   + "_view, ";
      } else {
        const char* arg_var_name = generated_all[&a_];
        this->var_name_ = name_gen.generate();
        generated[&a_] = this->var_name_.c_str();
        res.body = generate_body(row_index_name, col_index_name, view_handled,
                                 arg_var_name);
      }
    } else {
      this->var_name_ = generated[&a_];
    }
    return res;
  }

  /**
   * Generates kernel code for main body of this expression.
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @param view_handled whether whether caller already handled matrix view
   * @param arg_var_name variable name used for kernel arguments this is using
   * @return part of kernel with code for this expression
   */
  inline std::string generate_body(const std::string& row_index_name,
                                   const std::string& col_index_name,
                                   const bool view_handled,
                                   const char* arg_var_name) const {
    std::string type = type_str<Scalar>();
    if (view_handled) {
      return type + " " + var_name_ + " = " + arg_var_name + "_global["
             + row_index_name + " + " + arg_var_name + "_rows * "
             + col_index_name + "];\n";
    } else {
      return type + " " + var_name_ + " = 0;"
                 " if (!((!contains_nonzero(" + arg_var_name
                 + "_view, LOWER) && " + col_index_name + " < " + row_index_name
                 + ") || (!contains_nonzero(" + arg_var_name
                 + "_view, UPPER) && " + col_index_name + " > " + row_index_name
                 + "))) {" + var_name_ + " = " + arg_var_name + "_global["
                 + row_index_name + " + " + arg_var_name + "_rows * "
                 + col_index_name + "];}\n";
    }
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
    if (generated_all.count(&a_) == 0) {
      this->var_name_ = name_gen.generate();
    } else {
      this->var_name_ = generated_all[&a_];
    }
    kernel_parts res = generate_lhs(row_index_name, col_index_name);

    if (generated_all.count(&a_) == 0) {
      generated_all[&a_] = this->var_name_.c_str();
    } else {
      res.args = "";
    }
    return res;
  }

  /**
   * Generates kernel code for this expression if it appears on the left hand
   * side of an assignment.
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @return part of kernel with code for this expressions
   */
  inline kernel_parts generate_lhs(const std::string& row_index_name,
                                   const std::string& col_index_name) const {
    kernel_parts res;
    std::string type = type_str<Scalar>();
    res.args = "__global " + type + "* " + var_name_ + "_global, int "
               + var_name_ + "_rows, int " + var_name_ + "_view, ";
    res.body = var_name_ + "_global[" + row_index_name + " + " + var_name_
               + "_rows * " + col_index_name + "]";
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
    if (generated_all.count(&a_) == 0) {
      generated_all[&a_] = "";
      kernel.setArg(arg_num++, a_.buffer());
      kernel.setArg(arg_num++, a_.rows());
      kernel.setArg(arg_num++, a_.view());
    }
  }

  /**
   * Adds read event to the matrix used in this expression.
   * @param e the event to add
   */
  inline void add_read_event(cl::Event& e) const { a_.add_read_event(e); }

  /**
   * Adds write event to the matrix used in this expression.
   * @param e the event to add
   */
  inline void add_write_event(cl::Event& e) const { a_.add_write_event(e); }

  /**
   * Adds all read and write events on the matrix used by this expression to a
   * list and clears them from the matrix.
   * @param[out] events List of all events.
   */
  inline void get_clear_read_write_events(
      std::vector<cl::Event>& events) const {
    events.insert(events.end(), a_.read_events().begin(),
                  a_.read_events().end());
    events.insert(events.end(), a_.write_events().begin(),
                  a_.write_events().end());
    a_.clear_read_write_events();
  }

  /**
   * Adds all write events on the matrix used by this expression to a list.
   * @param[out] events List of all events.
   */
  inline void get_write_events(std::vector<cl::Event>& events) const {
    events.insert(events.end(), a_.write_events().begin(),
                  a_.write_events().end());
  }

  /**
   * Number of rows of a matrix that would be the result of evaluating this
   * expression.
   * @return number of rows
   */
  inline int rows() const { return a_.rows(); }

  /**
   * Number of columns of a matrix that would be the result of evaluating this
   * expression.
   * @return number of columns
   */
  inline int cols() const { return a_.cols(); }

  /**
   * View of a matrix that would be the result of evaluating this expression.
   * @return view
   */
  inline matrix_cl_view view() const { return a_.view(); }

  /**
   * Sets the view of the matrix depending on which of its parts are written to.
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
    if (bottom_zero_diagonal <= top_diagonal && bottom_diagonal < 0) {
      a_.view(either(a_.view(), matrix_cl_view::Lower));
    } else if (bottom_zero_diagonal <= 1 - a_.rows()) {
      a_.view(both(a_.view(), matrix_cl_view::Upper));
    }
    if (top_zero_diagonal >= bottom_diagonal && top_diagonal > 0) {
      a_.view(either(a_.view(), matrix_cl_view::Upper));
    } else if (top_zero_diagonal >= a_.cols() - 1) {
      a_.view(both(a_.view(), matrix_cl_view::Lower));
    }
  }

  /**
   * Determine indices of extreme sub- and superdiagonals written.
   * @return pair of indices - bottom and top diagonal
   */
  inline std::pair<int, int> extreme_diagonals() const {
    int bottom = contains_nonzero(a_.view(), matrix_cl_view::Lower)
                     ? -a_.rows() + 1
                     : 0;
    int top = contains_nonzero(a_.view(), matrix_cl_view::Upper) ? a_.cols() - 1
                                                                 : 0;
    return {bottom, top};
  }

  /**
   * Evaluates the expression. \c load_ returns a const reference to stored
   * matrix_cl.
   * @return Result of the expression.
   */
  const T& eval() const { return a_; }

  /**
   * If needed resizes underlying matrix to desired number of rows and cols.
   * @param rows desired number of rows
   * @param cols desired number of columns
   */
  inline void check_assign_dimensions(int rows, int cols) const {
    if (a_.rows() != rows || a_.cols() != cols) {
      a_ = matrix_cl<Scalar>(rows, cols);
    }
  }

  /**
   * Collects data that is needed beside types to uniqly identify a kernel
   * generator expression.
   * @param[out] uids ids of unique matrix accesses
   * @param[in,out] id_map map from memory addresses to unique ids
   * @param[in,out] next_id neqt unique id to use
   */
  inline void get_unique_matrix_accesses(
      std::vector<int>& uids, std::unordered_map<const void*, int>& id_map,
      int& next_id) const {
    if (id_map.count(&a_) == 0) {
      id_map[&a_] = next_id;
      uids.push_back(next_id);
      next_id++;
    } else {
      uids.push_back(id_map[&a_]);
    }
  }
};

/** @}*/
}  // namespace math
}  // namespace stan

#endif
#endif
