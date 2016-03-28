#include <stan/math/fwd/arr.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix,value_of_rec) {
  using stan::math::value_of_rec;
  using std::vector;
  using stan::math::fvar;

  vector<double> a_vals;

  for (size_t i = 0; i < 10; ++i)
    a_vals.push_back(i + 1);

  vector<double> b_vals;

  for (size_t i = 10; i < 15; ++i)
    b_vals.push_back(i + 1);

  {
    vector<fvar<double> > a;
    a = stan::math::to_fvar(a_vals);
    vector<fvar<double> > b;
    b = stan::math::to_fvar(b_vals);

    vector<double> d_a = value_of_rec(a);
    vector<double> d_b = value_of_rec(b);

    for (int i = 0; i < 5; ++i)
      EXPECT_FLOAT_EQ(b[i].val_, d_b[i]);

    for (int i = 0; i < 10; ++i)
      EXPECT_FLOAT_EQ(a[i].val_, d_a[i]);
  }

  {
    vector<fvar<double> > a_vals_fd = stan::math::to_fvar(a_vals);
    vector<fvar<double> > b_vals_fd = stan::math::to_fvar(b_vals);
    vector<fvar<double> > zeros_fd(10);
    std::fill(zeros_fd.begin(), zeros_fd.end(), 0);

    vector<fvar<fvar<double> > > a;
    a = stan::math::to_fvar(a_vals_fd, zeros_fd);
    vector<fvar<fvar<double> > > b;
    b = stan::math::to_fvar(b_vals_fd, zeros_fd);

    vector<double> d_a = value_of_rec(a);
    vector<double> d_b = value_of_rec(b);

    for (int i = 0; i < 5; ++i)
      EXPECT_FLOAT_EQ(b[i].val_.val_, d_b[i]);

    for (int i = 0; i < 10; ++i)
      EXPECT_FLOAT_EQ(a[i].val_.val_, d_a[i]);
  }

}
