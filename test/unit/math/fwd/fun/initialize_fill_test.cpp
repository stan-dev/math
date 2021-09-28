#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(AgradFwdMatrixinitialize_fill, fd) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::initialize_fill;
  using std::vector;

  fvar<double> x;
  fvar<double> y = 10;
  initialize_fill(x, y);
  EXPECT_FLOAT_EQ(10.0, x.val_);

  std::vector<fvar<double> > z(2);
  fvar<double> a = 15;
  initialize_fill(z, a);
  EXPECT_FLOAT_EQ(15.0, z[0].val_);
  EXPECT_FLOAT_EQ(15.0, z[1].val_);
  EXPECT_EQ(2U, z.size());

  Matrix<fvar<double>, Dynamic, Dynamic> m(2, 3);
  initialize_fill(m, fvar<double>(12));
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j)
      EXPECT_FLOAT_EQ(12.0, m(i, j).val_);

  Matrix<fvar<double>, Dynamic, 1> rv(3);
  initialize_fill(rv, fvar<double>(13));
  for (int i = 0; i < 3; ++i)
    EXPECT_FLOAT_EQ(13.0, rv(i).val_);

  Matrix<fvar<double>, 1, Dynamic> v(4);
  initialize_fill(v, fvar<double>(22));
  for (int i = 0; i < 4; ++i)
    EXPECT_FLOAT_EQ(22.0, v(i).val_);

  vector<vector<fvar<double> > > d(3, vector<fvar<double> >(2));
  initialize_fill(d, fvar<double>(54));
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 2; ++j)
      EXPECT_FLOAT_EQ(54, d[i][j].val_);
}
TEST(AgradFwdMatrixinitialize_fill, fd2) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::initialize_fill;
  Matrix<fvar<double>, Dynamic, 1> y(3);
  initialize_fill(y, 3.0);
  EXPECT_EQ(3, y.size());
  EXPECT_FLOAT_EQ(3.0, y[0].val_);
}

TEST(AgradFwdMatrixinitialize_fill, ffd) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::initialize_fill;
  using std::vector;

  fvar<fvar<double> > x;
  fvar<fvar<double> > y = 10;
  initialize_fill(x, y);
  EXPECT_FLOAT_EQ(10.0, x.val_.val_);

  std::vector<fvar<fvar<double> > > z(2);
  fvar<fvar<double> > a = 15;
  initialize_fill(z, a);
  EXPECT_FLOAT_EQ(15.0, z[0].val_.val_);
  EXPECT_FLOAT_EQ(15.0, z[1].val_.val_);
  EXPECT_EQ(2U, z.size());

  Matrix<fvar<fvar<double> >, Dynamic, Dynamic> m(2, 3);
  initialize_fill(m, fvar<fvar<double> >(12));
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j)
      EXPECT_FLOAT_EQ(12.0, m(i, j).val_.val_);

  Matrix<fvar<fvar<double> >, Dynamic, 1> rv(3);
  initialize_fill(rv, fvar<fvar<double> >(13));
  for (int i = 0; i < 3; ++i)
    EXPECT_FLOAT_EQ(13.0, rv(i).val_.val_);

  Matrix<fvar<fvar<double> >, 1, Dynamic> v(4);
  initialize_fill(v, fvar<fvar<double> >(22));
  for (int i = 0; i < 4; ++i)
    EXPECT_FLOAT_EQ(22.0, v(i).val_.val_);

  vector<vector<fvar<fvar<double> > > > d(3, vector<fvar<fvar<double> > >(2));
  initialize_fill(d, fvar<fvar<double> >(54));
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 2; ++j)
      EXPECT_FLOAT_EQ(54, d[i][j].val_.val_);
}
TEST(AgradFwdMatrixinitialize_fill, ffd2) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::initialize_fill;
  Matrix<fvar<fvar<double> >, Dynamic, 1> y(3);
  initialize_fill(y, 3.0);
  EXPECT_EQ(3, y.size());
  EXPECT_FLOAT_EQ(3.0, y[0].val_.val_);
}
