#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

TEST(AgradRevMatrix,mv_trace) {
  using stan::math::trace;
  using stan::math::matrix_v;

  matrix_v a(2,2);
  a << -1.0, 2.0, 
    5.0, 10.0;
  
  AVEC x = createAVEC(a(0,0), a(0,1), a(1,0), a(1,1));

  AVAR s = trace(a);
  EXPECT_FLOAT_EQ(9.0,s.val());
  
  VEC g = cgradvec(s,x);
  EXPECT_FLOAT_EQ(1.0, g[0]);
  EXPECT_FLOAT_EQ(0.0, g[1]);
  EXPECT_FLOAT_EQ(0.0, g[2]);
  EXPECT_FLOAT_EQ(1.0, g[3]);
}  

TEST(AgradRevMatrix, check_varis_on_stack) {
  stan::math::matrix_v a(2,2);
  a << -1.0, 2.0, 
    5.0, 10.0;
  test::check_varis_on_stack(stan::math::trace(a));
}
