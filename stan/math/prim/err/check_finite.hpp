#ifndef STAN_MATH_PRIM_ERR_CHECK_FINITE_HPP
#define STAN_MATH_PRIM_ERR_CHECK_FINITE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/is_scal_finite.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/err/throw_domain_error_vec.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {
template <typename T_y, bool is_vec>
struct finite {
  static void check(const char* function, const char* name, const T_y& y) {
    if (!is_scal_finite(y)) {
      throw_domain_error(function, name, y, "is ", ", but must be finite!");
    }
  }
};

template <typename T_y>
struct finite<T_y, true> {
  static void check(const char* function, const char* name, const T_y& y) {
    for (size_t n = 0; n < size(y); n++) {
      if (!is_scal_finite(stan::get(y, n))) {
        throw_domain_error_vec(function, name, y, n, "is ",
                               ", but must be finite!");
      }
    }
  }
};
}  // namespace internal

/**
 * Check if <code>y</code> is finite.
 * This function is vectorized and will check each element of
 * <code>y</code>.
 * @tparam T_y Type of y
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @throw <code>domain_error</code> if y is infinity, -infinity, or NaN
 */
template <typename T_y>
inline void check_finite(const char* function, const char* name, const T_y& y) {
  internal::finite<T_y, is_vector_like<T_y>::value>::check(function, name, y);
}

/*
 * Return <code>true</code> is the specified matrix is finite.
 *
 * @tparam T type of elements in the matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param function name of function (for error messages)
 * @param name variable name (for error messages)
 * @param y matrix to test
 * @return <code>true</code> if the matrix is finite
 **/
namespace internal {
template <typename T, int R, int C>
struct finite<Eigen::Matrix<T, R, C>, true> {
  static void check(const char* function, const char* name,
                    const Eigen::Matrix<T, R, C>& y) {
    if (!value_of(y).allFinite()) {
      for (int n = 0; n < y.size(); ++n) {
        if (!std::isfinite(value_of_rec(y(n)))) {
          throw_domain_error_vec(function, name, y, n, "is ",
                                 ", but must be finite!");
        }
      }
    }
  }
};

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
