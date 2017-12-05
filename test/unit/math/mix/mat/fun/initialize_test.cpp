#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>
#include <vector>

using stan::math::fvar;
using stan::math::var;

TEST(AgradMixMatrixInitialize, fv) {
  using stan::math::initialize;
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  fvar<var>  x;
  fvar<var>  y = 10;
  initialize(x, y);
  EXPECT_FLOAT_EQ(10.0, x.val_.val());

  std::vector<fvar<var> > z(2);
  fvar<var>  a = 15;
  initialize(z, a);
  EXPECT_FLOAT_EQ(15.0, z[0].val_.val());
  EXPECT_FLOAT_EQ(15.0, z[1].val_.val());
  EXPECT_EQ(2U, z.size());

  Matrix<fvar<var> , Dynamic, Dynamic> m(2, 3);
  initialize(m, fvar<var> (12));
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j)
      EXPECT_FLOAT_EQ(12.0, m(i, j).val_.val());

  Matrix<fvar<var> , Dynamic, 1> rv(3);
  initialize(rv, fvar<var> (13));
  for (int i = 0; i < 3; ++i)
    EXPECT_FLOAT_EQ(13.0, rv(i).val_.val());

  Matrix<fvar<var> , 1, Dynamic> v(4);
  initialize(v, fvar<var> (22));
  for (int i = 0; i < 4; ++i)
    EXPECT_FLOAT_EQ(22.0, v(i).val_.val());

  vector<vector<fvar<var> > > d(3, vector<fvar<var> >(2));
  initialize(d, fvar<var> (54));
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 2; ++j)
      EXPECT_FLOAT_EQ(54, d[i][j].val_.val());
}
TEST(AgradMixMatrixInitialize, fv2) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::initialize;
  Matrix<fvar<var> , Dynamic, 1> y(3);
  initialize(y, 3.0);
  EXPECT_EQ(3, y.size());
  EXPECT_FLOAT_EQ(3.0, y[0].val_.val());
}

TEST(AgradMixMatrixInitialize, ffv) {
  using stan::math::initialize;
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  fvar<fvar<var> >  x;
  fvar<fvar<var> >  y = 10;
  initialize(x, y);
  EXPECT_FLOAT_EQ(10.0, x.val_.val_.val());

  std::vector<fvar<fvar<var> > > z(2);
  fvar<fvar<var> >  a = 15;
  initialize(z, a);
  EXPECT_FLOAT_EQ(15.0, z[0].val_.val_.val());
  EXPECT_FLOAT_EQ(15.0, z[1].val_.val_.val());
  EXPECT_EQ(2U, z.size());

  Matrix<fvar<fvar<var> > , Dynamic, Dynamic> m(2, 3);
  initialize(m, fvar<fvar<var> > (12));
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j)
      EXPECT_FLOAT_EQ(12.0, m(i, j).val_.val_.val());

  Matrix<fvar<fvar<var> > , Dynamic, 1> rv(3);
  initialize(rv, fvar<fvar<var> > (13));
  for (int i = 0; i < 3; ++i)
    EXPECT_FLOAT_EQ(13.0, rv(i).val_.val_.val());

  Matrix<fvar<fvar<var> > , 1, Dynamic> v(4);
  initialize(v, fvar<fvar<var> > (22));
  for (int i = 0; i < 4; ++i)
    EXPECT_FLOAT_EQ(22.0, v(i).val_.val_.val());

  vector<vector<fvar<fvar<var> > > > d(3, vector<fvar<fvar<var> > >(2));
  initialize(d, fvar<fvar<var> > (54));
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 2; ++j)
      EXPECT_FLOAT_EQ(54, d[i][j].val_.val_.val());
}
TEST(AgradMixMatrixInitialize, ffv2) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::initialize;
  Matrix<fvar<fvar<var> > , Dynamic, 1> y(3);
  initialize(y, 3.0);
  EXPECT_EQ(3, y.size());
  EXPECT_FLOAT_EQ(3.0, y[0].val_.val_.val());
}
