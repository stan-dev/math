#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_BINARY_FOO_FUN_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_BINARY_FOO_FUN_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_binary.hpp>
#include <boost/math/tools/promotion.hpp>
#include <stan/math/prim/scal/err/check_less_or_equal.hpp>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * Example of how to define a functor for a vectorized function.
     * The example includes a constrained version of pow().
     */
    struct binary_foo_fun {
      template <typename T1, typename T2>
      static inline typename boost::math::tools::promote_args<T1, T2>::type
      fun(const T1& x, const T2& y) {
        using std::pow;
        stan::math::check_less_or_equal("binary_foo_fun vectorize",
                                        "x", x, 10);
        stan::math::check_less_or_equal("binary_foo_fun vectorize",
                                        "y", y, 10);
        return pow(x, y);
      }
    };

    template <typename T1, typename T2>
    inline typename apply_scalar_binary<binary_foo_fun, T1, T2>::return_t
    binary_foo(const T1& x, const T2& y) {
      return apply_scalar_binary<binary_foo_fun, T1, T2>::apply(x, y);
    }
  }
}

#endif
