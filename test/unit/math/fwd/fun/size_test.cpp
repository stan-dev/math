#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(AgradFwdMatrixSize, fvar_double) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::size;
  using std::vector;

  vector<fvar<double> > y(6);
  EXPECT_EQ(6, stan::math::size(y));

  vector<Matrix<fvar<double>, Dynamic, Dynamic> > z(7);
  EXPECT_EQ(7, stan::math::size(z));

  vector<Matrix<fvar<double>, Dynamic, 1> > a(8);
  EXPECT_EQ(8, stan::math::size(a));

  vector<Matrix<fvar<double>, 1, Dynamic> > b(9);
  EXPECT_EQ(9, stan::math::size(b));

  vector<vector<fvar<double> > > c(10);
  EXPECT_EQ(10, stan::math::size(c));

  vector<vector<fvar<double> > > ci(10);
  EXPECT_EQ(10, stan::math::size(ci));

  vector<vector<Matrix<fvar<double>, Dynamic, Dynamic> > > d(11);
  EXPECT_EQ(11, stan::math::size(d));

  vector<vector<Matrix<fvar<double>, 1, Dynamic> > > e(12);
  EXPECT_EQ(12, stan::math::size(e));

  vector<vector<Matrix<fvar<double>, Dynamic, 1> > > f(13);
  EXPECT_EQ(13, stan::math::size(f));

  vector<vector<vector<Matrix<fvar<double>, Dynamic, 1> > > > g(14);
  EXPECT_EQ(14, stan::math::size(g));
}

TEST(AgradFwdMatrixSize, fvar_fvar_double) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::size;
  using std::vector;

  vector<fvar<fvar<double> > > y(6);
  EXPECT_EQ(6, stan::math::size(y));

  vector<Matrix<fvar<fvar<double> >, Dynamic, Dynamic> > z(7);
  EXPECT_EQ(7, stan::math::size(z));

  vector<Matrix<fvar<fvar<double> >, Dynamic, 1> > a(8);
  EXPECT_EQ(8, stan::math::size(a));

  vector<Matrix<fvar<fvar<double> >, 1, Dynamic> > b(9);
  EXPECT_EQ(9, stan::math::size(b));

  vector<vector<fvar<fvar<double> > > > c(10);
  EXPECT_EQ(10, stan::math::size(c));

  vector<vector<fvar<fvar<double> > > > ci(10);
  EXPECT_EQ(10, stan::math::size(ci));

  vector<vector<Matrix<fvar<fvar<double> >, Dynamic, Dynamic> > > d(11);
  EXPECT_EQ(11, stan::math::size(d));

  vector<vector<Matrix<fvar<fvar<double> >, 1, Dynamic> > > e(12);
  EXPECT_EQ(12, stan::math::size(e));

  vector<vector<Matrix<fvar<fvar<double> >, Dynamic, 1> > > f(13);
  EXPECT_EQ(13, stan::math::size(f));

  vector<vector<vector<Matrix<fvar<fvar<double> >, Dynamic, 1> > > > g(14);
  EXPECT_EQ(14, stan::math::size(g));
}
