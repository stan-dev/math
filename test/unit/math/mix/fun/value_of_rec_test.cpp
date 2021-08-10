#include <stan/math/mix.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(AgradMix, value_of_rec) {
  using stan::math::fvar;
  using stan::math::value_of_rec;
  using stan::math::var;

  fvar<var> fv_a(5.0);
  fvar<fvar<var> > ffv_a(5.0);
  fvar<fvar<fvar<fvar<fvar<var> > > > > fffffv_a(5.0);

  EXPECT_FLOAT_EQ(5.0, value_of_rec(fv_a));
  EXPECT_FLOAT_EQ(5.0, value_of_rec(ffv_a));
  EXPECT_FLOAT_EQ(5.0, value_of_rec(fffffv_a));
}

TEST(MathMatrixMixArr, value_of_rec) {
  using stan::math::fvar;
  using stan::math::value_of_rec;
  using stan::math::var;
  using std::vector;

  vector<fvar<fvar<var> > > a;
  for (size_t i = 0; i < 10; ++i)
    a.push_back(fvar<fvar<var> >(i + 1));

  vector<fvar<fvar<var> > > b;
  for (size_t i = 10; i < 15; ++i)
    b.push_back(fvar<fvar<var> >(i + 1));

  vector<double> d_a = value_of_rec(a);
  vector<double> d_b = value_of_rec(b);

  for (int i = 0; i < 5; ++i)
    EXPECT_FLOAT_EQ(b[i].val_.val_.val(), d_b[i]);

  for (int i = 0; i < 10; ++i)
    EXPECT_FLOAT_EQ(a[i].val_.val_.val(), d_a[i]);
}

TEST(AgradMixMatrix, value_of_rec) {
  using stan::math::fvar;
  using stan::math::value_of_rec;
  using stan::math::var;
  using std::vector;

  vector<double> a_vals;

  for (size_t i = 0; i < 10; ++i)
    a_vals.push_back(i + 1);

  vector<double> b_vals;

  for (size_t i = 10; i < 15; ++i)
    b_vals.push_back(i + 1);

  Eigen::Matrix<fvar<var>, 2, 5> fv_a;
  ::fill(a_vals, fv_a);
  Eigen::Matrix<fvar<var>, 5, 1> fv_b;
  ::fill(b_vals, fv_b);

  Eigen::Matrix<fvar<fvar<var> >, 2, 5> ffv_a;
  ::fill(a_vals, ffv_a);
  Eigen::Matrix<fvar<fvar<var> >, 5, 1> ffv_b;
  ::fill(b_vals, ffv_b);

  Eigen::MatrixXd d_fv_a = value_of_rec(fv_a);
  Eigen::MatrixXd d_fv_b = value_of_rec(fv_b);
  Eigen::MatrixXd d_ffv_a = value_of_rec(ffv_a);
  Eigen::MatrixXd d_ffv_b = value_of_rec(ffv_b);

  for (Eigen::Index i = 0; i < 5; ++i) {
    EXPECT_FLOAT_EQ(b_vals[i], d_fv_b(i));
    EXPECT_FLOAT_EQ(b_vals[i], d_ffv_b(i));
  }

  for (Eigen::Index i = 0; i < 2; ++i)
    for (Eigen::Index j = 0; j < 5; ++j) {
      EXPECT_FLOAT_EQ(a_vals[j * 2 + i], d_fv_a(i, j));
      EXPECT_FLOAT_EQ(a_vals[j * 2 + i], d_ffv_a(i, j));
    }
}
