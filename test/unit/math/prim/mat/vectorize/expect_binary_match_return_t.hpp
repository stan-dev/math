#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_MATCH_RETURN_T_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_EXPECT_BINARY_MATCH_RETURN_T_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_binary.hpp>
#include <boost/type_traits/is_same.hpp>
#include <gtest/gtest.h>

namespace stan {

  namespace test {

    /**
     * Tests that the return type of a vectorized function with a
     * specified argument type is a specified type.
     *
     * @tparam F Functor defining function to test.
     * @tparam T_result_expected  Expected result type.
     * @tparam T_arg Argument type.
     */
    template <typename F, typename T_result_expected, 
              typename T_arg1, typename T_arg2>
    void expect_binary_match_return_t() {
      using stan::math::apply_scalar_binary;
      typedef typename apply_scalar_binary<F, T_arg1, T_arg2>::
      return_t result_t;
      EXPECT_TRUE((boost::is_same<T_result_expected, result_t>::value));
    }

  }
}
#endif


