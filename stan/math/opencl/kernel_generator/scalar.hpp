#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_CONSTANT_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_CONSTANT_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation.hpp>
#include <string>
#include <type_traits>
#include <set>

namespace stan {
namespace math {

/**
 * Represents a scalar in kernel generator expressions.
 * @tparam T type of the scalar
 */
template<typename T>
class scalar__ : public operation<scalar__<T>, T> {
public:
  static_assert(std::is_arithmetic<T>::value, "std::is_arithmetic<T> must be true for scalars!");
  using ReturnScalar = T;
  using base = operation<scalar__<T>, T>;
  using base::var_name;
  using base::instance;

  /**
   * Constructor
   * @param a scalar value
   */
  explicit scalar__(const T a) : a_(a) {}

  /**
 * generates kernel code for this and nested expressions.
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
      res.args = type_str<T>::name + " " + var_name + ", ";
      return res;
    }
    else {
      return {};
    }
  }

  /**
 * Sets kernel arguments for this and nested expressions.
 * @param[in,out] generated set of expressions that already set their kernel arguments
 * @param kernel kernel to set arguments on
 * @param[in,out] arg_num consecutive number of the first argument to set. This is incremented for each argument set by this function.
 */
  inline void set_args(std::set<int>& generated, cl::Kernel& kernel, int& arg_num) const {
    kernel.setArg(arg_num++, a_);
  }

  /**
 * Adds event for any matrices used by this expression.
 * @param e the event to add
 */
  inline void add_event(cl::Event& e) const {

  }

  /**
 * Number of rows of a matrix that would be the result of evaluating this expression.
 * @return number of rows
 */
  inline int rows() const {
    return base::dynamic;
  }

  /**
 * Number of columns of a matrix that would be the result of evaluating this expression.
 * @return number of columns
 */
  inline int cols() const {
    return base::dynamic;
  }

  /**
   * View of a matrix that would be the result of evaluating this expression.
   * @return view
   */
  inline matrix_cl_view view() const {
    return matrix_cl_view::Entire;
  }

private:
  T a_;
};

}
}

#endif
#endif
