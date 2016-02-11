#ifndef TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_MATCH_RETURN_T_HPP
#define TEST_UNIT_MATH_FWD_MAT_VECTORIZE_EXPECT_MATCH_RETURN_T_HPP

#include <stan/math/fwd/mat/vectorize/apply_scalar_unary.hpp>
#include <test/unit/math/prim/mat/vectorize/foo_fun.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <gtest/gtest.h>


namespace stan {
  namespace test {
    template <typename T1, typename T2>
    void expect_match_return_t() {
      using stan::math::foo_fun;
      using stan::math::apply_scalar_unary;
      using boost::is_convertible;

      typedef typename apply_scalar_unary<foo_fun, T2>::return_t result_t;
  
      bool T1_to_T2
        = is_convertible<T1, result_t>::value;
      EXPECT_TRUE(T1_to_T2);

      bool T1_from_T2
        = is_convertible<result_t, T1>::value;
      EXPECT_TRUE(T1_from_T2);
    }

  }
}
#endif
