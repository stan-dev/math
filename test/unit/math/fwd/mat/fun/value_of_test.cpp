#include <stan/math/fwd/mat.hpp>
#include <gtest/gtest.h>

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
