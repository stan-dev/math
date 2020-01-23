#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_LOAD_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_LOAD_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl_lhs.hpp>
#include <type_traits>
#include <string>
#include <utility>
#include <set>

namespace stan {
namespace math {

/**
 * Represents an access to a \c matrix_cl in kernel generator expressions
 * @tparam T \c matrix_cl
 */
template <typename T>
class load_
    : public operation_cl_lhs<load_<T>,
                              typename std::remove_reference_t<T>::type> {
 protected:
  T a_;

 public:
  using Scalar = typename std::remove_reference_t<T>::type;
  using base = operation_cl<load_<T>, Scalar>;
  using base::var_name;
  static_assert(std::is_base_of<matrix_cl<Scalar>,
                                typename std::remove_reference_t<T>>::value,
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
   * generates kernel code for this expression.
   * @param i row index variable name
   * @param j column index variable name
   * @return part of kernel with code for this expression
   */
  inline kernel_parts generate(const std::string& i,
                               const std::string& j) const {
    kernel_parts res{};
    std::string type = type_str<Scalar>();
    res.body = type + " " + var_name + " = 0;"
               " if (!((!contains_nonzero(" + var_name + "_view, LOWER) && "
               + j + " < " + i + ") || (!contains_nonzero(" + var_name +
               "_view, UPPER) && " + j + " > " + i + "))) {"
               + var_name + " = " + var_name + "_global[" + i + " + " +
               var_name + "_rows * " + j + "];}\n";
    res.args = "__global " + type + "* " + var_name + "_global, int " + var_name
               + "_rows, int " + var_name + "_view, ";
    return res;
  }

  /**
   * generates kernel code for this expression if it appears on the left hand
   * side of an assigment.
   * @param i row index variable name
   * @param j column index variable name
   * @return part of kernel with code for this expressions
   */
  inline kernel_parts generate_lhs(const std::string& i,
                                   const std::string& j) const {
    kernel_parts res;
    std::string type = type_str<Scalar>();
    res.args = "__global " + type + "* " + var_name + "_global, int " + var_name
               + "_rows, int " + var_name + "_view, ";
    res.body
        = var_name + "_global[" + i + " + " + var_name + "_rows * " + j + "]";
    return res;
  }

  /**
   * Sets kernel arguments for this expression.
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
   * Sets view of the matrix depending on which part is written.
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
    if (bottom_diagonal < 0) {
      a_.view(either(a_.view(), matrix_cl_view::Lower));
    } else if (bottom_zero_diagonal <= 1 - a_.rows()) {
      a_.view(both(a_.view(), matrix_cl_view::Upper));
    }
    if (top_diagonal > 0) {
      a_.view(either(a_.view(), matrix_cl_view::Upper));
    } else if (top_zero_diagonal >= a_.cols() - 1) {
      a_.view(both(a_.view(), matrix_cl_view::Lower));
    }
  }

  /**
   * Determine index of bottom diagonal written.
   * @return number of columns
   */
  inline int bottom_diagonal() const {
    return contains_nonzero(a_.view(), matrix_cl_view::Lower) ? -a_.rows() + 1
                                                              : 0;
  }

  /**
   * Determine index of top diagonal written.
   * @return number of columns
   */
  inline int top_diagonal() const {
    return contains_nonzero(a_.view(), matrix_cl_view::Upper) ? a_.cols() - 1
                                                              : 0;
  }

  /**
   * Evaluates the expression. \c load_ returns a const reference to stored
   * matrix_cl.
   * @return Result of the expression.
   */
  const T& eval() const { return a_; }
};

}  // namespace math
}  // namespace stan

#endif
#endif
