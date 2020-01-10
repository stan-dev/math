#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

TEST(MathPrimMat, vec_of_rvec_vec) {
  std::vector<Eigen::Matrix<double, -1, 1>> x(3);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(3, 1);
    x[i] << 1, 1, 1;
  }
  std::vector<double> vec(3);
  vec[0] = 1.0;
  vec[1] = 2.0;
  vec[2] = 3.0;
  std::vector<Eigen::Matrix<double, -1, 1>> out(3);
  out = stan::math::divide_columns(x, vec);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_FLOAT_EQ(out[i][j], 1.0 / (j + 1.0))
          << "index : ( " << i << ",  " << j << ")";
    }
  }
}

TEST(MathPrimMat, err_size_match) {
  std::vector<Eigen::Matrix<double, -1, 1>> x(3);
  for (size_t i = 0; i < x.size(); ++i) {
    x[i].resize(3, 1);
    x[i] << 1 * i, 2 * i, 3 * i;
  }
  std::vector<double> vec(2);
  vec[0] = 2.0;
  vec[1] = 3.0;
  EXPECT_THROW(stan::math::divide_columns(x, vec), std::invalid_argument);
}
