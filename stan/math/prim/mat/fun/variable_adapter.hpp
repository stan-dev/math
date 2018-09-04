#ifndef STAN_MATH_PRIM_MAT_FUN_VARIABLE_ADAPTER_HPP
#define STAN_MATH_PRIM_MAT_FUN_VARIABLE_ADAPTER_HPP

#include <stan/math/prim/scal/functor/apply.hpp>
#include <stan/math/prim/scal/err/check_less.hpp>
#include <vector>
#include <cstddef>
#include <iostream>

namespace stan {
namespace math {

/**
 * variable_adapter exports a flat, single-index interface to variables of
 * a given type in a variadic list of arguments provided in the constructor
 *
 * variable_adapter behaves like a 1d data structure with length equal to
 * the number of scalars of type T in the list of input arguments
 *
 * @tparam T Type of element to index
 * @tparam Targs Types of input arguments to index over
 */
template <typename T, typename... Targs>
class variable_adapter {
 protected:
  std::tuple<Targs...> args_;

  size_t size_;

  /**
   * count_T_impl is the implementation of count_T
   */
  template <typename... Pargs, int RowType, int ColType>
  size_t count_T_impl(size_t count, const Eigen::Matrix<T, RowType, ColType>& x,
                      const Pargs&... args) {
    return count_T_impl(count + x.size(), args...);
  }

  template <typename... Pargs>
  size_t count_T_impl(size_t count, const std::vector<T>& x,
                      const Pargs&... args) {
    return count_T_impl(count + x.size(), args...);
  }

  template <typename... Pargs>
  size_t count_T_impl(size_t count, const T& x, const Pargs&... args) {
    return count_T_impl(count + 1, args...);
  }

  template <typename R, typename... Pargs>
  size_t count_T_impl(size_t count, const R& x, const Pargs&... args) {
    return count_T_impl(count, args...);
  }

  size_t count_T_impl(size_t count) { return count; }

  /**
   * count_T counts the number of scalars of type T in the input argument list
   *
   * @tparam Pargs Types of input arguments
   * @return Number of scalars of type T in input
   */
  template <typename... Pargs>
  size_t count_T(const Pargs&... args) {
    return count_T_impl(0, args...);
  }

  /**
   * Tail of get recursion. This indicates an out of bounds error
   *
   * @return always throws
   * @throw domain_error always
   */
  T& get(size_t i) {
    std::ostringstream message;
    message << "variable_adapter::get: i (the index) points " << i
            << " places beyond the end of this variable_adapter";
    throw std::domain_error(message.str());
  }

  /**
   * get the ith element of arg. If i is greater than arg.size(), get
   * the (i - arg.size())th element of the remaining args
   *
   * @tparam RowType Eigen row type of arg
   * @tparam ColType Eigen column type of arg
   * @tparam Pargs Types of the rest of the input arguments to process
   * @return Reference to ith T in args_
   */
  template <int RowType, int ColType, typename... Pargs>
  T& get(size_t i, Eigen::Matrix<T, RowType, ColType>& arg, Pargs&... args) {
    if (i < static_cast<size_t>(arg.size()))
      return arg(i);
    else
      return get(i - arg.size(), args...);
  }

  /**
   * get the ith element of arg. If i is greater than arg.size(), get
   * the (i - arg.size())th element of the remaining args
   *
   * @tparam Pargs Types of the rest of the input arguments to process
   * @return Reference to ith T in args_
   */
  template <typename... Pargs>
  T& get(size_t i, std::vector<T>& arg, Pargs&... args) {
    if (i < arg.size())
      return arg[i];
    else
      return get(i - arg.size(), args...);
  }

  /**
   * If i == 0, arg is the element we're looking for, otherwise, get
   * the (i - 1)th element of the remaining args
   *
   * @tparam Pargs Types of the rest of the input arguments to process
   * @return Reference to ith T in args_
   */
  template <typename... Pargs>
  T& get(size_t i, T& arg, Pargs&... args) {
    if (i == 0)
      return arg;
    else
      return get(i - 1, args...);
  }

  /**
   * No T here, recursively call get on the remaining arguments
   *
   * @tparam R Type of arg (is not indexed)
   * @tparam Pargs Types of the rest of the input arguments to process
   * @return Reference to ith T in args_
   */
  template <typename R, typename... Pargs>
  auto& get(size_t i, R& arg, Pargs&... args) {
    return get(i, args...);
  }

 public:
  /**
   * Construct a variable_adapter from the given list of inputs. The inputs
   * are copied by value into a tuple member variable args_, so writing to
   * variable_adapter does not modify the original variables.
   *
   * @param args Arguments to adapt to a 1D interface
   */
  variable_adapter(const Targs&... args)
      : args_(std::make_tuple(args...)), size_(count_T(args...)) {}

  /**
   * Get a tuple of args of the current values. These are returned by
   * value, not reference
   *
   * @return Tuple of args
   */
  std::tuple<Targs...> get_args() const { return args_; }

  /**
   * Return the number of Ts that were passed to the constructor
   *
   * @return Number of indexable elements
   */
  size_t size() const { return size_; }

  /**
   * Index the variable_adapter. Elements are returned by reference, so can
   * be modified
   *
   * @param i index into adapter
   * @return Reference to the ith T in args_
   */
  auto& operator()(size_t i) {
    check_less("variable_adapter::operator()", "i", i, size_);

    return apply(
        [ i, this ](auto&... args) -> auto& { return this->get(i, args...); },
        args_);
  }
};

/**
 * make_variable_adapter makes it easy to construct a variable_adapter
 * variable (template arguments are determined automatically)
 *
 * @tparam T Type of variable to index in args
 * @tparam Targs Types of input arguments
 * @param args Input arguments to build a 1d interface around
 * @return Newly constructed variable_adapter
 */
template <typename T, typename... Targs>
auto make_variable_adapter(const Targs&... args) {
  return variable_adapter<T, Targs...>(args...);
}

}  // namespace math
}  // namespace stan
#endif
