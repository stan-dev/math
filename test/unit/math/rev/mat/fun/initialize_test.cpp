#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix, initializeVar) {
  using stan::math::var;
  using stan::math::initialize;
  var a;
  var b = 10;
  initialize(a, b);           // template 1
  EXPECT_FLOAT_EQ(10, a.val());

  initialize(a, 5);          // template 2
  EXPECT_FLOAT_EQ(5, a.val());

  initialize(a, 13.2);      // template 2
  EXPECT_FLOAT_EQ(13.2, a.val());
}

TEST(MathMatrix, initMatrix) {
  using stan::math::var;
  using stan::math::initialize;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<var, Dynamic, Dynamic> mvar(2, 3);
  initialize(mvar, 2.3);     // template 3, 1
  for (int i = 0; i < mvar.size(); ++i)
    EXPECT_FLOAT_EQ(mvar(i).val(), 2.3);
}
