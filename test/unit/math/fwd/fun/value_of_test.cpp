#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(AgradFwd, value_of) {
  using stan::math::fvar;
  using stan::math::value_of;

  fvar<double> a = 5.0;
  EXPECT_FLOAT_EQ(5.0, value_of(a));
  // make sure all work together
  EXPECT_FLOAT_EQ(5.0, value_of(5.0));
  EXPECT_FLOAT_EQ(5.0, value_of(5));
}

TEST(AgradFwd, value_of_nan) {
  using stan::math::fvar;
  using stan::math::value_of;
  double nan = std::numeric_limits<double>::quiet_NaN();

  fvar<double> a = nan;
  EXPECT_TRUE(stan::math::is_nan(value_of(a)));
  EXPECT_TRUE(stan::math::is_nan(value_of(nan)));
}

TEST(MathMatrixFwdArr, value_of) {
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

TEST(AgradFwdMatrix, value_of) {
  using stan::math::fvar;
  using stan::math::value_of;
  using std::vector;

  Eigen::Matrix<double, 2, 5> a_vals;
  for (int i = 0; i < 10; ++i)
    a_vals(i) = i + 1;

  Eigen::Matrix<double, 5, 1> b_vals;
  for (int i = 0; i < 5; ++i)
    b_vals(i) = 10 + i;

  Eigen::Matrix<double, 2, 5> a = a_vals;
  Eigen::Matrix<double, 5, 1> b = b_vals;

  Eigen::Matrix<fvar<double>, 2, 5> fd_a;
  fd_a = stan::math::to_fvar(a_vals);
  Eigen::Matrix<fvar<double>, 5, 1> fd_b;
  fd_b = stan::math::to_fvar(b_vals);

  Eigen::Matrix<fvar<double>, 2, 5> fd_a_zeros;
  Eigen::Matrix<fvar<double>, 5, 1> fd_b_zeros;
  Eigen::Matrix<fvar<fvar<double> >, 2, 5> ffd_a;
  ffd_a = stan::math::to_fvar(fd_a, fd_a_zeros);
  Eigen::Matrix<fvar<fvar<double> >, 5, 1> ffd_b;
  ffd_b = stan::math::to_fvar(fd_b, fd_b_zeros);

  Eigen::MatrixXd d_fd_a = value_of(fd_a);
  Eigen::MatrixXd d_fd_b = value_of(fd_b);
  Eigen::Matrix<fvar<double>, -1, -1> d_ffd_a = value_of(ffd_a);
  Eigen::Matrix<fvar<double>, -1, -1> d_ffd_b = value_of(ffd_b);

  for (int i = 0; i < 5; ++i) {
    EXPECT_FLOAT_EQ(b(i), d_fd_b(i));
    EXPECT_FLOAT_EQ(b(i), d_ffd_b(i).val_);
  }

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 5; ++j) {
      EXPECT_FLOAT_EQ(a(i, j), d_fd_a(i, j));
      EXPECT_FLOAT_EQ(a(i, j), d_ffd_a(i, j).val_);
    }
}
