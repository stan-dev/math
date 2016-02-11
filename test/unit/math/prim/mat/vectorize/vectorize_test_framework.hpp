#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_VECTORIZE_TEST_FRAMEWORK_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_VECTORIZE_TEST_FRAMEWORK_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <test/unit/math/prim/mat/vectorize/foo_fun.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <gtest/gtest.h>


namespace stan {
  namespace test {

    template <double val>
    struct test_vectorize {

      static inline void test_function {
  
        using stan::math::foo;

        double exp_val = foo(val);
        EXPECT_FLOAT_EQ(std::exp(val), exp_val);

        EXPECT_FLOAT_EQ(std::exp(val), foo(val));
      }
    };

    template <typename T>
    struct test_vectorize {

      static inline void test_function {

        using stan::math::foo;
         
        
      } 
    };
  }
}
#endif
