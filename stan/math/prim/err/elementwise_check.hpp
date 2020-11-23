#ifndef STAN_MATH_PRIM_ERR_ELEMENTWISE_ERROR_CHECKER_HPP
#define STAN_MATH_PRIM_ERR_ELEMENTWISE_ERROR_CHECKER_HPP

#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/fun/get.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <string>
#include <sstream>
#include <vector>

namespace stan {
namespace math {
namespace internal {

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

/**
 * Pipes given arguments into a stringstream. Integral arguments (indices) are
 * increased by 1 before.
 *
 * @tparam Arg0 type of the first argument
 * @tparam Args types of remaining arguments
 * @param ss stringstream to pipe arguments in
 * @param arg0 the first argument
 * @param args remining arguments
 */
inline void pipe_in(std::stringstream& ss) {}
template <typename Arg0, typename... Args>
inline void pipe_in(std::stringstream& ss, Arg0 arg0, const Args... args) {
  if (std::is_integral<Arg0>::value) {
    ss << arg0 + 1;
  } else {
    ss << arg0;
  }
  pipe_in(ss, args...);
}

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
template <typename F, typename T, typename... Indexings,
          require_stan_scalar_t<T>* = nullptr>
inline void elementwise_check(const F& is_good, const char* function,
                              const char* name, const T& x, const char* must_be,
                              const Indexings... indexings) {
  if (unlikely(!is_good(x))) {
    std::stringstream ss{};
    ss << function << ": " << name;
    internal::pipe_in(ss, indexings...);
    ss << "is " << x << ", but must be " << must_be << "!";
    throw std::domain_error(ss.str());
  }
}
template <typename F, typename T, typename... Indexings,
          require_eigen_t<T>* = nullptr,
          std::enable_if_t<static_cast<bool>(Eigen::internal::traits<T>::Flags&(
              Eigen::LinearAccessBit | Eigen::DirectAccessBit))>* = nullptr>
inline void elementwise_check(const F& is_good, const char* function,
                              const char* name, const T& x, const char* must_be,
                              const Indexings... indexings) {
  for (size_t i = 0; i < x.size(); i++) {
    auto scal = value_of(x.coeff(i));
    if (unlikely(!is_good(scal))) {
      std::stringstream ss{};
      ss << function << ": " << name;
      internal::pipe_in(ss, indexings...);
      if (is_eigen_vector<T>::value) {
        ss << "[" << i + 1 << "] is " << scal << ", but must be " << must_be
           << "!";
      } else if (Eigen::internal::traits<T>::Flags & Eigen::RowMajorBit) {
        ss << "[" << i / x.rows() + 1 << ", " << i % x.rows() + 1 << "] is "
           << scal << ", but must be " << must_be << "!";
      } else {
        ss << "[" << i % x.cols() + 1 << ", " << i / x.cols() + 1 << "] is "
           << scal << ", but must be " << must_be << "!";
      }
      throw std::domain_error(ss.str());
    }
  }
}
template <
    typename F, typename T, typename... Indexings,
    require_eigen_t<T>* = nullptr,
    std::enable_if_t<!(Eigen::internal::traits<T>::Flags
                       & (Eigen::LinearAccessBit | Eigen::DirectAccessBit))
                     && !(Eigen::internal::traits<T>::Flags
                          & Eigen::RowMajorBit)>* = nullptr>
inline void elementwise_check(const F& is_good, const char* function,
                              const char* name, const T& x, const char* must_be,
                              const Indexings... indexings) {
  for (size_t i = 0; i < x.rows(); i++) {
    for (size_t j = 0; j < x.cols(); j++) {
      auto scal = value_of(x.coeff(i, j));
      if (unlikely(!is_good(scal))) {
        std::stringstream ss{};
        ss << function << ": " << name;
        internal::pipe_in(ss, indexings...);
        ss << "[" << i + 1 << ", " << j + 1 << "] is " << scal
           << ", but must be " << must_be << "!";
        throw std::domain_error(ss.str());
      }
    }
  }
}
template <
    typename F, typename T, typename... Indexings,
    require_eigen_t<T>* = nullptr,
    std::enable_if_t<!(Eigen::internal::traits<T>::Flags
                       & (Eigen::LinearAccessBit | Eigen::DirectAccessBit))
                     && static_cast<bool>(Eigen::internal::traits<T>::Flags
                                          & Eigen::RowMajorBit)>* = nullptr>
inline void elementwise_check(const F& is_good, const char* function,
                              const char* name, const T& x, const char* must_be,
                              const Indexings... indexings) {
  for (size_t j = 0; j < x.cols(); j++) {
    for (size_t i = 0; i < x.rows(); i++) {
      auto scal = value_of(x.coeff(i, j));
      if (unlikely(!is_good(scal))) {
        std::stringstream ss{};
        ss << function << ": " << name;
        internal::pipe_in(ss, indexings...);
        ss << "[" << i + 1 << ", " << j + 1 << "] is " << scal
           << ", but must be " << must_be << "!";
        throw std::domain_error(ss.str());
      }
    }
  }
}
template <typename F, typename T, typename... Indexings,
          require_std_vector_t<T>* = nullptr>
inline void elementwise_check(const F& is_good, const char* function,
                              const char* name, const T& x, const char* must_be,
                              const Indexings... indexings) {
  for (size_t j = 0; j < x.size(); j++) {
    elementwise_check(is_good, function, name, x[j], must_be, indexings..., "[",
                      j, "]");
  }
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
