#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_FOO_FUN_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_FOO_FUN_HPP

#include <cmath>

namespace stan {

  namespace math {
    
    // mock class static function definition for exp()
    struct foo_fun {
      template <typename T>
      static inline T fun(const T& x) {
        using std::exp;
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
