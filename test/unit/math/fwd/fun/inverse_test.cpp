#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(AgradFwdMatrixDims, matrix_fd) {
  using stan::math::matrix_fd;
  using stan::math::inverse;

  matrix_fd t1 = Eigen::MatrixXd::Random(3, 3);
  matrix_fd res = inverse(t1);
}
