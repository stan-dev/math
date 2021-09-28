#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <vector>

TEST(AgradRevMatrix, fill) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fill;
  using std::vector;

  stan::math::var x;
  stan::math::var y = 10;
  fill(x, y);
  EXPECT_FLOAT_EQ(10.0, x.val());

  std::vector<stan::math::var> z(2);
  stan::math::var a = 15;
  fill(z, a);
  EXPECT_FLOAT_EQ(15.0, z[0].val());
  EXPECT_FLOAT_EQ(15.0, z[1].val());
  EXPECT_EQ(2U, z.size());

  Matrix<stan::math::var, Dynamic, Dynamic> m(2, 3);
  fill(m, stan::math::var(12));
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j)
      EXPECT_FLOAT_EQ(12.0, m(i, j).val());

  Matrix<stan::math::var, Dynamic, 1> rv(3);
  fill(rv, stan::math::var(13));
  for (int i = 0; i < 3; ++i)
    EXPECT_FLOAT_EQ(13.0, rv(i).val());

  Matrix<stan::math::var, 1, Dynamic> v(4);
  fill(v, stan::math::var(22));
  for (int i = 0; i < 4; ++i)
    EXPECT_FLOAT_EQ(22.0, v(i).val());

  vector<vector<stan::math::var>> d(3, vector<stan::math::var>(2));
  fill(d, stan::math::var(54));
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 2; ++j)
      EXPECT_FLOAT_EQ(54, d[i][j].val());
}
TEST(AgradRevMatrix, fillDouble) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fill;
  Matrix<double, Dynamic, 1> y(3);
  fill(y, 3.0);
  EXPECT_EQ(3, y.size());
  EXPECT_FLOAT_EQ(3.0, y[0]);
}

TEST(AgradRevMatrix, fillVarMatDouble) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fill;
  using stan::math::sum;
  using stan::math::var_value;
  Matrix<double, Dynamic, 1> y_val(3);
  var_value<Matrix<double, Dynamic, 1>> y(y_val);
  fill(y, 3.0);
  EXPECT_EQ(3, y.size());
  EXPECT_FLOAT_EQ(3.0, y.val()[0]);
  sum(y).grad();
  for (Eigen::Index i = 0; i < y.size(); ++i) {
    EXPECT_FLOAT_EQ(y.val()(i), y_val(i));
    EXPECT_FLOAT_EQ(y.adj()(i), 0);
  }
}

TEST(AgradRevMatrix, fillVarMatVar) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fill;
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  Matrix<double, Dynamic, 1> y_val(3);
  var_value<Matrix<double, Dynamic, 1>> y(y_val);
  var z(3.0);
  fill(y, z);
  EXPECT_EQ(3, y.size());
  EXPECT_FLOAT_EQ(3.0, y.val()[0]);
  sum(y).grad();
  EXPECT_FLOAT_EQ(z.val(), 3.0);
  EXPECT_FLOAT_EQ(z.adj(), 3.0);
  for (Eigen::Index i = 0; i < y.size(); ++i) {
    EXPECT_FLOAT_EQ(y.val()(i), y_val(i));
    EXPECT_FLOAT_EQ(y.adj()(i), 0);
  }
}
