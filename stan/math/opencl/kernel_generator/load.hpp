#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_LOAD_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_LOAD_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation.hpp>
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
template<typename T>
class load__ : public operation<load__<T>, typename std::remove_reference_t<T>::type> {
public:
  using ReturnScalar = typename std::remove_reference_t<T>::type;
  using base = operation<load__<T>, ReturnScalar>;
  using base::var_name;
  using base::instance;
  static_assert(std::is_base_of<matrix_cl < ReturnScalar>, typename std::remove_reference_t<T>>::value, "load__: argument a must be a matrix_cl<T>!");
  static_assert(std::is_arithmetic<ReturnScalar>::value, "load__: T in \"matrix_cl<T> a\" argument must be an arithmetic type!");

  /**
   * Constructor
   * @param a \c matrix_cl
   */
  explicit load__(T&& a) : a_(std::forward<T>(a)) {}

  /**
 * generates kernel code for this expression.
 * @param ng name generator for this kernel
 * @param[in,out] generated set of already generated operations
 * @param i row index variable name
 * @param j column index variable name
 * @return part of kernel with code for this and nested expressions
 */
  inline kernel_parts generate(name_generator& ng, std::set<int>& generated, const std::string& i, const std::string& j) const {
    if (generated.count(instance) == 0) {
      generated.insert(instance);
      var_name = ng.generate();
      kernel_parts res;
      std::string type = type_str<ReturnScalar>::name;
      res.body = type + " " + var_name + " = 0;"
                                         " if (!((!contains_nonzero(" + var_name + "_view, LOWER) && " + j + " < " + i + ") || (!contains_nonzero(" + var_name + "_view, UPPER) && " + j + " > " + i +
                 "))) {"
                 + var_name + " = " + var_name + "_global[" + i + " + " + var_name + "_rows * " + j + "];}\n";
      res.args = "__global " + type + "* " + var_name + "_global, int " + var_name + "_rows, int " + var_name + "_view, ";
      return res;
    }
    else {
      return {};
    }
  }

  /**
   * generates kernel code for this expression if it appears on the left hand side of an assigment.
   * @param ng name generator for this kernel
   * @param[in,out] generated set of already generated operations
   * @param i row index variable name
   * @param j column index variable name
   * @return part of kernel with code for this expressions
   */
  inline kernel_parts generate_lhs(name_generator& ng, std::set<int>& generated, const std::string& i, const std::string& j) const {
    kernel_parts res;
    if (generated.count(instance) == 0) {
      generated.insert(instance);
      var_name = ng.generate();
      std::string type = type_str<ReturnScalar>::name;
      res.args = "__global " + type + "* " + var_name + "_global, int " + var_name + "_rows, int " + var_name + "_view, ";
    }
    res.body = var_name + "_global[" + i + " + " + var_name + "_rows * " + j + "]";
    return res;
  }

  /**
 * Sets kernel arguments for this expression.
 * @param[in,out] generated set of expressions that already set their kernel arguments
 * @param kernel kernel to set arguments on
 * @param[in,out] arg_num consecutive number of the first argument to set. This is incremented for each argument set by this function.
 */
  inline void set_args(std::set<int>& generated, cl::Kernel& kernel, int& arg_num) const {
    if (generated.count(instance) == 0) {
      generated.insert(instance);
      kernel.setArg(arg_num++, a_.buffer());
      kernel.setArg(arg_num++, a_.rows());
      kernel.setArg(arg_num++, a_.view());
    }
  }

  /**
   * Adds event to the matrix used in this expression.
   * @param e the event to add
   */
  inline void add_event(cl::Event& e) const {
    a_.add_read_event(e);
  }

  /**
   * Adds write event to the matrix used in this expression.
   * @param e the event to add
   */
  inline void add_write_event(cl::Event& e) const {
    a_.add_write_event(e);
  }

  /**
   * Number of rows of a matrix that would be the result of evaluating this expression.
   * @return number of rows
   */
  inline int rows() const {
    return a_.rows();
  }

  /**
 * Number of columns of a matrix that would be the result of evaluating this expression.
 * @return number of columns
 */
  inline int cols() const {
    return a_.cols();
  }

  /**
   * View of a matrix that would be the result of evaluating this expression.
   * @return view
   */
  inline matrix_cl_view view() const {
    return a_.view();
  }

protected:
  T a_;
};

}
}

#endif
#endif
