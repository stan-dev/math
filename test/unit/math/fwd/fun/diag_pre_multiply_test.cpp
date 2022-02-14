#include <stan/math/fwd.hpp>
#include <test/unit/util.hpp>

#include <gtest/gtest.h>

TEST(MathFunRev, dpm_ug) {
  using stan::math::matrix_fd;
  using stan::math::vector_fd;
  using stan::math::diag_pre_multiply;

  Eigen::MatrixXd dm1 = Eigen::MatrixXd::Random(3,3);
  matrix_fd vm1 = dm1;

  Eigen::VectorXd dv2 = Eigen::VectorXd::Random(3);
  vector_fd vv2 = dv2;

  auto out = diag_pre_multiply(dv2, vm1);
  out = diag_pre_multiply(vv2, vm1);
  out = diag_pre_multiply(vv2, dm1);
}

TEST(MathFunFwd, dpm_ug_ffd) {
  using stan::math::matrix_ffd;
  using stan::math::vector_ffd;
  using stan::math::diag_pre_multiply;

  Eigen::MatrixXd dm1 = Eigen::MatrixXd::Random(3,3);
  matrix_ffd vm1 = dm1;

  Eigen::VectorXd dv2 = Eigen::VectorXd::Random(3);
  vector_ffd vv2 = dv2;

  auto out = diag_pre_multiply(dv2, vm1);
  out = diag_pre_multiply(vv2, vm1);
  out = diag_pre_multiply(vv2, dm1);
}