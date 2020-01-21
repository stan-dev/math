#ifndef STAN_MATH_PRIM_ERR_ELEMENTWISE_ERROR_CHECKER_HPP
#define STAN_MATH_PRIM_ERR_ELEMENTWISE_ERROR_CHECKER_HPP

#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/fun/get.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <string>
#include <sstream>
#include <vector>

namespace stan {
namespace math {
namespace internal {

template <typename T_is_good, typename T_exception>
struct Checker {
  const T_is_good& is_good;
  const char* function;
  const char* name;
  const char* suffix;

  void raise_error_impl(std::stringstream& ss) { throw T_exception{ss.str()}; }

  template <typename... M>
  void raise_error_impl(std::stringstream& ss, const char* message,
                        const M&... messages) {
    ss << message;
    raise_error_impl(ss, messages...);
  }

  template <typename... M>
  void raise_error_impl(std::stringstream& ss, double value_of_x,
                        const M&... messages) {
    ss << value_of_x;
    raise_error_impl(ss, messages...);
  }

  template <typename... M>
  void raise_error_impl(std::stringstream& ss, size_t index,
                        const M&... messages) {
    ss << index + 1;
    raise_error_impl(ss, messages...);
  }

  template <typename... M>
  void raise_error(const M&... messages) {
    std::stringstream ss{};
    raise_error_impl(ss, messages...);
  }

  template <typename T, typename = require_stan_scalar_t<T>,
            typename... T_indices>
  void check(const T& x, T_indices... indices) {
    double xd = value_of_rec(x);
    if (!is_good(xd))
      raise_error(function, ": ", name, indices..., " is ", xd, suffix);
  }

  template <typename T, typename = require_eigen_vector_t<T>, typename = void,
            typename... T_indices>
  void check(const T& x, T_indices... indices) {
    for (size_t i = 0; i < stan::math::size(x); ++i)
      check(x(i), indices..., "[", i, "]");
  }

  template <typename T, typename... T_indices>
  void check(const std::vector<T>& x, T_indices... indices) {
    for (size_t i = 0; i < stan::math::size(x); ++i)
      check(x[i], indices..., "[", i, "]");
  }

  template <typename T, typename... T_indices>
  void check(const Eigen::DenseBase<T>& x, T_indices... indices) {
    for (size_t n = 0; n < x.cols(); ++n)
      for (size_t m = 0; m < x.rows(); ++m)
        check(x(m, n), indices..., "[row=", m, ", col=", n, "]");
  }
};

template <typename T_is_good>
struct Iser {
  const T_is_good& is_good;

  template <typename T, typename = require_stan_scalar_t<T>>
  bool is(const T& x) {
    return is_good(value_of_rec(x));
  }

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
 * Check that the predicate holds for the value of x, working elementwise on
 * containers. If x is a container, check each element inside x, recursively.
 * @tparam T_is_good type of is_good predicate
 * @tparam T_x type of x
 * @param is_good predicate to check, must accept doubles and produce bools
 * @param function function name (for error messages)
 * @param name variable name (for error messages)
 * @param x variable to check, can be a scalar, a container of scalars, a
 * container of containers of scalars, etc
 * @param suffix message to print at end of error message
 * @throw <code>std::domain_error</code> if is_good returns false for the value
 * of any element in x the predicate
 */
template <typename T_is_good, typename T_x>
void elementwise_check(const T_is_good& is_good, const char* function,
                       const char* name, const T_x& x, const char* suffix) {
  internal::Checker<T_is_good, std::domain_error>{is_good, function, name,
                                                  suffix}
      .check(x);
}

/**
 * Like elementwise_check, but indicate failure by returning false instead of
 * by throwing.
 */
template <typename T_is_good, typename T_x>
bool elementwise_is(const T_is_good& is_good, const T_x& x) {
  return internal::Iser<T_is_good>{is_good}.is(x);
}

}  // namespace math
}  // namespace stan
#endif
