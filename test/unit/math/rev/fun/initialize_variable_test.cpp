#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <vector>

TEST(AgradRevMatrix, initializeVariable) {
  using stan::math::initialize_variable;
  using std::vector;

  using Eigen::Dynamic;
  using Eigen::Matrix;

  stan::math::var a;
  initialize_variable(a, stan::math::var(1.0));
  EXPECT_FLOAT_EQ(1.0, a.val());

  std::vector<stan::math::var> b(3);
  initialize_variable(b, stan::math::var(2.0));
  EXPECT_EQ(3U, b.size());
  EXPECT_FLOAT_EQ(2.0, b[0].val());
  EXPECT_FLOAT_EQ(2.0, b[1].val());
  EXPECT_FLOAT_EQ(2.0, b[2].val());

  vector<std::vector<stan::math::var>> c(4, std::vector<stan::math::var>(3));
  initialize_variable(c, stan::math::var(3.0));
  for (size_t m = 0; m < c.size(); ++m)
    for (size_t n = 0; n < c[0].size(); ++n)
      EXPECT_FLOAT_EQ(3.0, c[m][n].val());

  Matrix<stan::math::var, Dynamic, Dynamic> aa(5, 7);
  initialize_variable(aa, stan::math::var(4.0));
  for (int m = 0; m < aa.rows(); ++m)
    for (int n = 0; n < aa.cols(); ++n)
      EXPECT_FLOAT_EQ(4.0, aa(m, n).val());

  Matrix<stan::math::var, Dynamic, 1> bb(5);
  initialize_variable(bb, stan::math::var(5.0));
  for (int m = 0; m < bb.size(); ++m)
    EXPECT_FLOAT_EQ(5.0, bb(m).val());

  Matrix<stan::math::var, 1, Dynamic> cc(12);
  initialize_variable(cc, stan::math::var(7.0));
  for (int m = 0; m < cc.size(); ++m)
    EXPECT_FLOAT_EQ(7.0, cc(m).val());

  Matrix<stan::math::var, Dynamic, Dynamic> init_val(3, 4);
  vector<Matrix<stan::math::var, Dynamic, Dynamic>> dd(5, init_val);
  initialize_variable(dd, stan::math::var(11.0));
  for (size_t i = 0; i < dd.size(); ++i)
    for (int m = 0; m < dd[0].rows(); ++m)
      for (int n = 0; n < dd[0].cols(); ++n)
        EXPECT_FLOAT_EQ(11.0, dd[i](m, n).val());
}
