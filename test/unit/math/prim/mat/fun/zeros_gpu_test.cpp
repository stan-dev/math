#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix,zero_m_exception_pass) {
  stan::math::matrix_gpu m(1,1);
  
  EXPECT_NO_THROW(stan::math::identity(m));
  EXPECT_NO_THROW(stan::math::zeros(m, stan::math::UPPER));
  EXPECT_NO_THROW(stan::math::zeros(m, stan::math::LOWER));

  // TODO(Rok/Steve): Don't know if we should throw an error here?
  //stan::math::matrix_gpu m0(0,0);
  //EXPECT_NO_THROW(stan::math::identity(m0));
  
}
