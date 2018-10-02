#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_FOO_FUN_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_FOO_FUN_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/err/check_less_or_equal.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

/**
 * Example of how to define a functor for a vectorized function.
 * The example includes a constrained, punctured version of exp().
 */
struct foo_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::exp;
    stan::math::check_less_or_equal("foo_fun vectorize", "x", x, 5);
    if (x == 0)
      return std::numeric_limits<double>::quiet_NaN();
    return exp(x);
  }
};

template <typename T>
inline typename apply_scalar_unary<foo_fun, T>::return_t foo(const T& x) {
  return apply_scalar_unary<foo_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
