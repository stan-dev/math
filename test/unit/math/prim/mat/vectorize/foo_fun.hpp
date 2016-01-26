#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_FOO_FUN_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_FOO_FUN_HPP

#include <cmath>
#include <stan/math/prim/scal/err/check_less_or_equal.hpp>
#include <stan/math.hpp>

namespace stan {

  namespace math {
    
    // mock class static function definition for exp()
    struct foo_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using std::exp;
        stan::math::check_less_or_equal("vectorize", "x", x, 5);
        return exp(x);
      }
    };

    template <typename T>
    inline typename apply_scalar_unary<foo_fun, T>::return_t
    foo(const T& x) {
      return apply_scalar_unary<foo_fun, T>::apply(x);
    }

  }
}

#endif
