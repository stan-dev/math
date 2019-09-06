#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_BINARY_OPERATOR_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_BINARY_OPERATOR_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/err/check_matching_dims.hpp>
#include <stan/math/opencl/multiply.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation.hpp>
#include <stan/math/opencl/kernel_generator/scalar.hpp>
#include <stan/math/opencl/kernel_generator/as_operation.hpp>
#include <stan/math/opencl/kernel_generator/is_valid_expression.hpp>
#include <string>
#include <type_traits>
#include <set>
#include <utility>

namespace stan {
namespace math {

/**
 * Represents a binary operation in kernel generator expressions.
 * @tparam Derived derived type
 * @tparam T_a type of first argument
 * @tparam T_b type of second argument
 */
template<typename Derived, typename T_a, typename T_b>
class binary_operation : public operation<Derived, typename std::common_type<typename std::remove_reference_t<T_a>::ReturnScalar, typename std::remove_reference_t<T_b>::ReturnScalar>::type> {
public:
  static_assert(std::is_base_of<operation_base,std::remove_reference_t<T_a>>::value, "binary_operation: a must be an operation!");
  static_assert(std::is_base_of<operation_base,std::remove_reference_t<T_b>>::value, "binary_operation: b must be an operation!");

  using ReturnScalar = typename std::common_type<typename std::remove_reference_t<T_a>::ReturnScalar, typename std::remove_reference_t<T_b>::ReturnScalar>::type;
  using base = operation<Derived, ReturnScalar>;
  using base::var_name;
  using base::instance;

  /**
   * Constructor
   * @param a first argument
   * @param b sedond argument
   * @param op operation
   */
  binary_operation(T_a&& a, T_b&& b, const std::string& op) : a_(std::forward<T_a>(a)), b_(std::forward<T_b>(b)), op_(op) {
      const std::string function = "binary_operator" + op;
      if(a.rows()!=base::dynamic && b.rows()!=base::dynamic) {
        check_size_match(function.c_str(), "Rows of ", "a", a.rows(), "rows of ", "b", b.rows());
      }
      if(a.cols()!=base::dynamic && b.cols()!=base::dynamic) {
        check_size_match(function.c_str(), "Columns of ", "a", a.cols(), "columns of ", "b", b.cols());
      }
  }

  /**
   * generates kernel code for this and nested expressions.
   * @param ng name generator for this kernel
   * @param[in,out] generated set of already generated operations
   * @param i row index variable name
   * @param j column index variable name
   * @return part of kernel with code for this and nested expressions
   */
  inline kernel_parts generate(name_generator& ng, std::set<int>& generated, const std::string& i, const std::string& j) const{
    if(generated.count(instance)==0) {
      kernel_parts a_parts = a_.generate(ng, generated, i, j);
      kernel_parts b_parts = b_.generate(ng, generated, i, j);
      generated.insert(instance);
      var_name = ng.generate();
      kernel_parts res;
      res.body = a_parts.body + b_parts.body + type_str<ReturnScalar>::name + " " + var_name + " = " + a_.var_name + op_ + b_.var_name + ";\n";
      res.args = a_parts.args + b_parts.args;
      return res;
    }
    else{
      return {};
    }
  }

  /**
   * Sets kernel arguments for this and nested expressions.
   * @param[in,out] generated set of expressions that already set their kernel arguments
   * @param kernel kernel to set arguments on
   * @param[in,out] arg_num consecutive number of the first argument to set. This is incremented for each argument set by this function.
   */
  inline void set_args(std::set<int>& generated, cl::Kernel& kernel, int& arg_num) const{
    if(generated.count(instance)==0) {
      generated.insert(instance);
      a_.set_args(generated, kernel, arg_num);
      b_.set_args(generated, kernel, arg_num);
    }
  }

  /**
   * Adds event for any matrices used by this or nested expressions.
   * @param e the event to add
   */
  inline void add_event(cl::Event& e) const{
    a_.add_event(e);
    b_.add_event(e);
  }

  /**
   * Number of rows of a matrix that would be the result of evaluating this expression.
   * @return number of rows
   */
  inline int rows() const{
    int a_rows = a_.rows();
    return a_rows == base::dynamic ? b_.rows() : a_rows;
  }

  /**
   * Number of columns of a matrix that would be the result of evaluating this expression.
   * @return number of columns
   */
  inline int cols() const{
    int a_cols = a_.cols();
    return a_cols == base::dynamic ? b_.cols() : a_cols;
  }

  /**
   * View of a matrix that would be the result of evaluating this expression.
   * @return view
   */
  inline matrix_cl_view view() const{
    return either(a_.view(), b_.view());
  }

protected:
  T_a a_;
  T_b b_;
  std::string op_;
};

/**
 * Represents addition in kernel generator expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 */
template<typename T_a, typename T_b>
class addition__ : public binary_operation<addition__<T_a, T_b>, T_a, T_b> {
public:
  /**
   * Constructor.
   * @param a first expression
   * @param b second expression
   */
  addition__(T_a&& a, T_b&& b) : binary_operation<addition__<T_a, T_b>, T_a, T_b>(std::forward<T_a>(a),std::forward<T_b>(b),"+") {}
};

/**
 * Addition of two kernel generator expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first argument
 * @param b second argument
 * @return Addition of given expressions
 */
template<typename T_a, typename T_b, typename = enable_if_all_valid_expressions<T_a, T_b>>
inline addition__<as_operation_t<T_a>,as_operation_t<T_b>> operator+(T_a&& a, T_b&& b) {
  return {as_operation(std::forward<T_a>(a)), as_operation(std::forward<T_b>(b))};
}

/**
 * Represents subtraction in kernel generator expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 */
template<typename T_a, typename T_b>
class subtraction__ : public binary_operation<subtraction__<T_a, T_b>, T_a, T_b> {
public:
  /**
 * Constructor.
 * @param a first expression
 * @param b second expression
 */
  subtraction__(T_a&& a, T_b&& b) : binary_operation<subtraction__<T_a, T_b>, T_a, T_b>(std::forward<T_a>(a),std::forward<T_b>(b),"-") {}
};

/**
 * Subtraction of two kernel generator expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 * @return Subtraction of given expressions
 */
template<typename T_a, typename T_b, typename = enable_if_all_valid_expressions<T_a, T_b>>
inline subtraction__<as_operation_t<T_a>,as_operation_t<T_b>> operator-(T_a&& a, T_b&& b){
  return {as_operation(std::forward<T_a>(a)), as_operation(std::forward<T_b>(b))};
}

/**
 * Represents element-wise multiplication in kernel generator expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 */
template<typename T_a, typename T_b, typename = enable_if_all_valid_expressions<T_a, T_b>>
class elewise_multiplication__ : public binary_operation<elewise_multiplication__<T_a, T_b>, T_a, T_b> {
public:
  /**
  * Constructor.
  * @param a first expression
  * @param b second expression
  */
  elewise_multiplication__(T_a&& a, T_b&& b) : binary_operation<elewise_multiplication__<T_a, T_b>, T_a, T_b>(std::forward<T_a>(a),std::forward<T_b>(b),"*") {}

  /**
 * View of a matrix that would be the result of evaluating this expression.
 * @return view
 */
  inline matrix_cl_view view()  const{
    using base = binary_operation<elewise_multiplication__<T_a, T_b>, T_a, T_b>;
    return both(base::a_.view(), base::b_.view());
  }
};

/**
 * Element-wise multiplication of two kernel generator expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 * @return Element-wise multiplication of given expressions
 */
template<typename T_a, typename T_b>
inline elewise_multiplication__<as_operation_t<T_a>,as_operation_t<T_b>> elewise_multiplication(T_a&& a, T_b&& b) {
  return {as_operation(std::forward<T_a>(a)), as_operation(std::forward<T_b>(b))};
}

/**
 * Multiplication of a scalar and a kernel generator expression.
 * @tparam T_a type of scalar
 * @tparam T_b type of expression
 * @param a scalar
 * @param b expression
 * @return Multiplication of given arguments
 */
template<typename T_a, typename T_b, typename = enable_if_arithmetic<T_a>, typename = enable_if_all_valid_expressions<T_b>>
inline elewise_multiplication__<scalar__<T_a>,as_operation_t<T_b>> operator*(T_a&& a, T_b&& b) {
  return {as_operation(std::forward<T_a>(a)), as_operation(std::forward<T_b>(b))};
}

/**
 * Multiplication of a kernel generator expression and a scalar.
 * @tparam T_a type of expression
 * @tparam T_b type of scalar
 * @param a expression
 * @param b scalar
 * @return Multiplication of given arguments
 */
template<typename T_a, typename T_b, typename = enable_if_all_valid_expressions<T_a>, typename = enable_if_arithmetic<T_b>>
inline elewise_multiplication__<as_operation_t<T_a>,scalar__<T_b>> operator*(T_a&& a, const T_b b) {
  return {as_operation(std::forward<T_a>(a)), as_operation(b)};
}

/**
 * Matrix multiplication of two kernel generator expressions. Evaluates both expressions before calculating the matrix product.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 * @return Matrix product of given arguments
 */
template<typename T_a, typename T_b, typename = enable_if_all_valid_expressions_and_none_scalar<T_a, T_b>>
inline matrix_cl<double> operator*(const T_a& a, const T_b& b){
  //no need for perfect forwarding as operations are evaluated
  return as_operation(a).eval() * as_operation(b).eval();
}

/**
 * Represents element-wise division in kernel generator expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 */
template<typename T_a, typename T_b>
class elewise_division__ : public binary_operation<elewise_division__<T_a, T_b>, T_a, T_b> {
public:
  /**
  * Constructor.
  * @param a first expression
  * @param b second expression
  */
  elewise_division__(T_a&& a, T_b&& b) : binary_operation<elewise_division__<T_a, T_b>, T_a, T_b>(std::forward<T_a>(a),std::forward<T_b>(b),"/") {}

  /**
 * View of a matrix that would be the result of evaluating this expression.
 * @return view
 */
  inline matrix_cl_view view()  const{
    using base = binary_operation<elewise_division__<T_a, T_b>, T_a, T_b>;
    return either(base::a_.view(), invert(base::b_.view()));
  }
};

/**
 * Element-wise division of two kernel generator expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 * @return Element-wise division of given expressions
 */
template<typename T_a, typename T_b, typename = enable_if_all_valid_expressions<T_a, T_b>>
inline elewise_division__<as_operation_t<T_a>,as_operation_t<T_b>> elewise_division(T_a&& a, T_b&& b) {
  return {as_operation(std::forward<T_a>(a)), as_operation(std::forward<T_b>(b))};
}

}
}
#endif
#endif
