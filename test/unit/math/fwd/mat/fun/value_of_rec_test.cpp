#include <stan/math/fwd/mat.hpp>
#include <gtest/gtest.h>

TEST(AgradMatrix, value_of_rec) {
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
