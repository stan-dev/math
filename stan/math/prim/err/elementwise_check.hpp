#ifndef STAN_MATH_PRIM_ERR_ELEMENTWISE_ERROR_CHECKER_HPP
#define STAN_MATH_PRIM_ERR_ELEMENTWISE_ERROR_CHECKER_HPP

#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/fun/get.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <string>
#include <sstream>
#include <vector>

namespace stan {
namespace math {
namespace internal {

/** Apply an error check to a container, signal failure by throwing.
 * Apply a predicate like is_positive to the double underlying every scalar in a
 * container, throw an exception if the predicate fails for any double.
 * @tparam F type of predicate
 * @tparam E type of exception thrown
 */
template <typename F, typename E>
class Checker {
  const F& is_good;
  const char* function;
  const char* name;
  const char* suffix;

  /**
   * Throw an exception of type `E`.
   * The error message is the string inside the provided stringstream.
   * @param ss stringstream containing error message
   * @throws `E`
   */
  void raise_error_ss(std::stringstream& ss) { throw E{ss.str()}; }

  /**
   * Throw an exception of type `E`.
   * The error message is the concatenation of the string inside the provided
   * stringstream with all the provided messages.
   * @tparam M types of first message
   * @tparam Ms types of other messages
   * @param ss stringstream to accumulate error message in.
   * @param message a message to append to `ss`
   * @param messages more messages to append
   * @throws `E`
   */
  template <typename M, typename... Ms>
  void raise_error_ss(std::stringstream& ss, const M& message,
                      const Ms&... messages) {
    ss << message;
    raise_error_ss(ss, messages...);
  }

  /**
   * Throw an exception of type `E`.
   * The error message is the concatenation of the provided messages.
   * @tparam Ms types of messages
   * @param messages a list of messages
   * @throws `E`
   */
  template <typename... Ms>
  void raise_error(const Ms&... messages) {
    std::stringstream ss{};
    raise_error_ss(ss, messages...);
  }

 public:
  /**
   * @param is_good predicate to check, must accept doubles and produce bools
   * @param function function name (for error messages)
   * @param name variable name (for error messages)
   * @param suffix message to print at end of error message
   */
  Checker(const F& is_good, const char* function, const char* name,
          const char* suffix)
      : is_good(is_good), function(function), name(name), suffix(suffix) {}

  /**
   * Check the scalar.
   * @tparam T type of scalar
   * @tparam Ms types of messages
   * @param x scalar
   * @param messages a list of messages to append to the error message
   * @throws `E` if  the scalar fails the error check
   */
  template <typename T, require_stan_scalar_t<T>* = nullptr, typename... Ms>
  void check(const T& x, Ms... messages) {
    double xd = value_of_rec(x);
    if (!is_good(xd))
      raise_error(function, ": ", name, messages..., " is ", xd, suffix);
  }

  /**
   * Check all the scalars inside the standard vector.
   * @tparam T type of vector
   * @tparam Ms types of messages
   * @param x vector
   * @param messages a list of messages to append to the error message
   * @throws `E` if any of the scalars fail the error check
   */
  template <typename T, require_std_vector_t<T>* = nullptr, typename... Ms>
  void check(const T& x, Ms... messages) {
    for (size_t i = 0; i < stan::math::size(x); ++i)
      check(x[i], messages..., "[", i + 1, "]");
  }

  /**
   * Check all the scalars inside an eigen vector.
   * @tparam T type of vector
   * @tparam Ms types of messages
   * @param x vector
   * @param messages a list of messages to append to the error message
   * @throws `E` if any of the scalars fail the error check
   */
  template <typename T, require_eigen_vector_t<T>* = nullptr, typename... Ms>
  void check(const T& x, Ms... messages) {
    for (size_t i = 0; i < stan::math::size(x); ++i)
      check(x.coeff(i), messages..., "[", i + 1, "]");
  }

  /**
   * Check all the scalars inside the `var_value<Matrix>`.
   * @tparam T type of vector
   * @tparam Ms types of messages
   * @param x vector
   * @param messages a list of messages to append to the error message
   * @throws `E` if any of the scalars fail the error check
   */
  template <typename T, require_var_matrix_t<T>* = nullptr, typename... Ms>
  void check(const T& x, Ms... messages) {
    check(x.val(), messages...);
  }

  /**
   * Check all the scalars inside the matrix.
   * @tparam Derived type of matrix
   * @tparam Ms types of messages
   * @param x matrix
   * @param messages a list of messages to append to the error message
   * @throws `E` if any of the scalars fail the error check
   */
  template <typename EigMat, require_eigen_matrix_dynamic_t<EigMat>* = nullptr,
            typename... Ms>
  void check(const EigMat& x, Ms... messages) {
    for (size_t n = 0; n < x.cols(); ++n)
      for (size_t m = 0; m < x.rows(); ++m)
        check(x.coeff(m, n), messages..., "[row=", m + 1, ", col=", n + 1, "]");
  }
};  // namespace internal

/** Apply an error check to a container, signal failure with `false`.
 * Apply a predicate like is_positive to the double underlying every scalar in a
 * container, producing true if the predicate holds everywhere and `false` if it
 * fails anywhere.
 * @tparam F type of predicate
 */
template <typename F>
class Iser {
  const F& is_good;

 public:
  /**
   * @param is_good predicate to check, must accept doubles and produce bools
   */
  explicit Iser(const F& is_good) : is_good(is_good) {}

  /**
   * Check the scalar.
   * @tparam T type of scalar
   * @param x scalar
   * @return `false` if the scalar fails the error check
   */
  template <typename T, typename = require_stan_scalar_t<T>>
  bool is(const T& x) {
    return is_good(value_of_rec(x));
  }

  /**
   * Check all the scalars inside the container.
   * @tparam T type of scalar
   * @param x container
   * @return `false` if any of the scalars fail the error check
   */
  template <typename T, typename = require_not_stan_scalar_t<T>,
            typename = void>
  bool is(const T& x) {
    for (size_t i = 0; i < stan::math::size(x); ++i)
      if (!is(stan::get(x, i)))
        return false;
    return true;
  }
};

}  // namespace internal

/**
 * Check that the predicate holds for the value of `x`, working elementwise on
 * containers. If `x` is a scalar, check the double underlying the scalar. If
 * `x` is a container, check each element inside `x`, recursively.
 * @tparam F type of predicate
 * @tparam T type of `x`
 * @param is_good predicate to check, must accept doubles and produce bools
 * @param function function name (for error messages)
 * @param name variable name (for error messages)
 * @param x variable to check, can be a scalar, a container of scalars, a
 * container of containers of scalars, etc
 * @param suffix message to print at end of error message
 * @throws `std::domain_error` if `is_good` returns `false` for the value
 * of any element in `x`
 */
template <typename F, typename T>
inline void elementwise_check(const F& is_good, const char* function,
                              const char* name, const T& x,
                              const char* suffix) {
  internal::Checker<F, std::domain_error>{is_good, function, name, suffix}
      .check(x);
}

/**
 * Check that the predicate holds for the value of `x`, working elementwise on
 * containers. If `x` is a scalar, check the double underlying the scalar. If
 * `x` is a container, check each element inside `x`, recursively.
 * @tparam F type of predicate
 * @tparam T type of `x`
 * @param is_good predicate to check, must accept doubles and produce bools
 * @param x variable to check, can be a scalar, a container of scalars, a
 * container of containers of scalars, etc
 * @return `false` if any of the scalars fail the error check
 */
template <typename F, typename T>
inline bool elementwise_is(const F& is_good, const T& x) {
  return internal::Iser<F>{is_good}.is(x);
}

}  // namespace math
}  // namespace stan
#endif
