#include <stan/math/rev/core.hpp>
#include <gtest/gtest.h>
#include <sstream>

TEST(AgradRev, insertion_operator) {
  stan::math::vari v(5);
  std::stringstream ss;
  ss << &v;
  EXPECT_EQ("5:0", ss.str());
}

TEST(AgradRev, long_double_test) {
  stan::math::vari_value<long double> v(5);
  std::stringstream ss;
  ss << &v;
  EXPECT_EQ("5:0", ss.str());
}

TEST(AgradRev, dense_matrix_vari) {
  using stan::math::vari_value;
  stan::math::vari_value<Eigen::MatrixXd> A(Eigen::MatrixXd::Random(3, 3));
}
