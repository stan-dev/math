#include <stan/math/rev/core.hpp>
#include <gtest/gtest.h>
#include <sstream>

TEST(AgradRev, insertion_operator) {
  stan::math::vari v(5);
  std::stringstream ss;
  ss << &v;
  EXPECT_EQ("5:0", ss.str());
}

TEST(AgradRev, eigen_obj) {
  using eigen_mat = Eigen::Matrix<double, -1, -1>;
  stan::math::vari_value<eigen_mat> A(eigen_mat::Random(10, 10));
  stan::math::vari_value<eigen_mat> B(eigen_mat::Random(10, 10), true);
  stan::math::vari_value<eigen_mat> C(eigen_mat::Random(10, 10), false);
  using eigen_arr = Eigen::Array<double, -1, -1>;
  stan::math::vari_value<eigen_arr> A_arr(eigen_arr::Random(10, 10));
  stan::math::vari_value<eigen_arr> B_arr(eigen_arr::Random(10, 10), true);
  stan::math::vari_value<eigen_arr> C_arr(eigen_arr::Random(10, 10), false);
  stan::math::recover_memory();
}
