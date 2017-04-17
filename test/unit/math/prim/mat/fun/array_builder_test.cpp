#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

using std::vector;
using Eigen::Dynamic;
using Eigen::Matrix;
using stan::math::array_builder;
using stan::math::to_matrix;
using stan::math::to_row_vector;

TEST(MathMatrix,arrayBuilder) {

  vector<Matrix<double, 1, Dynamic> > mEmpty
    = array_builder<Matrix<double, 1, Dynamic> >().array();
  EXPECT_EQ(0U,mEmpty.size());

  vector<double> v2 =
    array_builder<double>()
    .add(1)
    .add(2)
    .array();
  Matrix<double, 1, Dynamic> rv2 = to_row_vector(v2);
  vector<Matrix<double, 1, Dynamic> > v3rv2
    = array_builder<Matrix<double, 1, Dynamic> >().
    add(rv2).
    add(rv2).
    add(rv2).
    array();
  Matrix<double, 3, 2> m32 = to_matrix(v3rv2);

  EXPECT_EQ(6U,m32.size());
  EXPECT_FLOAT_EQ(1.0, m32(0,0));
  EXPECT_FLOAT_EQ(2.0, m32(0,1));
  EXPECT_FLOAT_EQ(1.0, m32(1,0));
  EXPECT_FLOAT_EQ(2.0, m32(1,1));
  EXPECT_FLOAT_EQ(1.0, m32(2,0));
  EXPECT_FLOAT_EQ(2.0, m32(2,1));

}

