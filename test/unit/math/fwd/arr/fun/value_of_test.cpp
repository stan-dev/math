#include <stan/math/fwd/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMatrix, value_of) {
  using stan::math::fvar;
  using stan::math::value_of;
  using std::vector;

  vector<double> a_vals;
  for (size_t i = 0; i < 10; ++i)
    a_vals.push_back(i + 1);

  vector<double> b_vals;
  for (size_t i = 10; i < 15; ++i)
    b_vals.push_back(i + 1);

  vector<fvar<double> > a_fd;
  a_fd = stan::math::to_fvar(a_vals);
  vector<fvar<double> > b_fd;
  b_fd = stan::math::to_fvar(b_vals);

  vector<double> d_a = value_of(a_fd);
  vector<double> d_b = value_of(b_fd);

  for (int i = 0; i < 5; ++i)
    EXPECT_FLOAT_EQ(b_fd[i].val_, d_b[i]);

  for (int i = 0; i < 10; ++i)
    EXPECT_FLOAT_EQ(a_fd[i].val_, d_a[i]);

  vector<fvar<double> > a_zeros(10);
  vector<fvar<fvar<double> > > a_ffd;
  a_ffd = stan::math::to_fvar(a_fd, a_zeros);

  vector<fvar<double> > b_zeros(5);
  vector<fvar<fvar<double> > > b_ffd;
  b_ffd = stan::math::to_fvar(b_fd, b_zeros);

  vector<fvar<double> > d_a_fd = value_of(a_ffd);
  vector<fvar<double> > d_b_fd = value_of(b_ffd);

  for (int i = 0; i < 5; ++i)
    EXPECT_FLOAT_EQ(b_ffd[i].val_.val_, d_b_fd[i].val_);

  for (int i = 0; i < 10; ++i)
    EXPECT_FLOAT_EQ(a_ffd[i].val_.val_, d_a_fd[i].val_);
}
