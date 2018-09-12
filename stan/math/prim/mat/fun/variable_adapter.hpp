#ifndef STAN_MATH_PRIM_MAT_FUN_VARIABLE_ADAPTER_HPP
#define STAN_MATH_PRIM_MAT_FUN_VARIABLE_ADAPTER_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/functor/apply.hpp>
#include <stan/math/prim/scal/err/check_less.hpp>
#include <cstddef>
#include <iostream>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

/**
 * Export a one-dimensional interface to scalar variables of a given
 * type in a variadic list of arguments provided in the constructor
 *
 * @tparam T Type of element to index
 * @tparam Targs Types of input arguments to index over
 */
template <typename T, typename... Targs>
class variable_adapter {
 protected:
  std::tuple<Targs...> args_;

  size_t size_;

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
   * Count the number of scalars of type T in the input argument list
   *
   * @tparam Pargs Types of input arguments
   * @return Number of scalars of type T in input
   */
  template <typename... Pargs>
  size_t count_T(const Pargs&... args) {
    return count_T_impl(0, args...);
  }

  /**
   * Tail of get recursion. This function always indicates an out of bounds
   * error
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
   * Return the ith element of arg if i is less than arg.size().
   * Otherwise return the (i - arg.size())th element of the
   * remaining args
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
   * Return the ith element of arg if i is less than arg.size().
   * Otherwise return the (i - arg.size())th element of the
   * remaining args
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
   * Return arg if i == 0. Otherwise, return
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
   * Return the ith element of the remaining arguments
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
   * The inputs are copied by value into a tuple member variable
   * args_, so writing to variable_adapter does not modify the
   * original variables.
   *
   * @param args Arguments to adapt to a 1D interface
   */
  explicit variable_adapter(const Targs&... args)
      : args_(std::make_tuple(args...)), size_(count_T(args...)) {}

  /**
   * Get a copy of the args_ member variable
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
   * Reference the elements of the variable_adapter
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
 * Construct a variable_adapter targetting type T from a variadic list of
 * arguments
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
