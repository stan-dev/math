#include <stan/math/fwd.hpp>
#include <stan/math/fwd/core.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(AgradFwd, value_of_rec) {
  using stan::math::fvar;

  fvar<double> a = 5.0;
  fvar<fvar<double> > ff_a(5.0);
  fvar<fvar<fvar<fvar<double> > > > ffff_a(5.0);
  EXPECT_FLOAT_EQ(5.0, value_of_rec(a));
  EXPECT_FLOAT_EQ(5.0, value_of_rec(ff_a));
  EXPECT_FLOAT_EQ(5.0, value_of_rec(ffff_a));
}

TEST(MathMatrixFwdArr, value_of_rec) {
  using stan::math::fvar;
  using stan::math::value_of_rec;
  using std::vector;

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

TEST(AgradMatrixFwd, value_of_rec) {
  using stan::math::fvar;
  using stan::math::value_of_rec;

  Eigen::Matrix<double, 2, 5> a;
  for (int i = 0; i < 10; ++i)
    a(i) = i + 1;

  Eigen::Matrix<double, 5, 1> b;
  for (int i = 0; i < 5; ++i)
    b(i) = 10 + i + 1;

  Eigen::Matrix<fvar<double>, 2, 5> fd_a;
  fd_a = stan::math::to_fvar(a);
  Eigen::Matrix<fvar<double>, 5, 1> fd_b;
  fd_b = stan::math::to_fvar(b);
  Eigen::Matrix<fvar<double>, 2, 5> fd_a_zeros;
  Eigen::Matrix<fvar<double>, 5, 1> fd_b_zeros;

  Eigen::Matrix<fvar<fvar<double> >, 2, 5> ffd_a;
  ffd_a = stan::math::to_fvar(fd_a, fd_a_zeros);
  Eigen::Matrix<fvar<fvar<double> >, 5, 1> ffd_b;
  ffd_b = stan::math::to_fvar(fd_b, fd_b_zeros);

  Eigen::MatrixXd d_fd_a = value_of_rec(fd_a);
  Eigen::MatrixXd d_fd_b = value_of_rec(fd_b);
  Eigen::MatrixXd d_ffd_a = value_of_rec(ffd_a);
  Eigen::MatrixXd d_ffd_b = value_of_rec(ffd_b);

  for (int i = 0; i < 5; ++i) {
    EXPECT_FLOAT_EQ(b(i), d_fd_b(i));
    EXPECT_FLOAT_EQ(b(i), d_ffd_b(i));
  }

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 5; ++j) {
      EXPECT_FLOAT_EQ(a(i, j), d_fd_a(i, j));
      EXPECT_FLOAT_EQ(a(i, j), d_ffd_a(i, j));
    }
}
