#include <stan/math/rev/core.hpp>
#include <test/unit/util.hpp>
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

TEST(AgradRev, dense_matrix_vari_block) {
  using stan::math::vari_value;
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(3, 3);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(3, 3);

  vari_value<Eigen::MatrixXd> A(a);
  EXPECT_MATRIX_EQ(a.block(0, 1, 2, 2), A.block(0, 1, 2, 2).val_);
  vari_value<Eigen::MatrixXd> B(a, b);
  EXPECT_MATRIX_EQ(a.block(0, 1, 2, 2), B.block(0, 1, 2, 2).val_);
  EXPECT_MATRIX_EQ(b.block(0, 1, 2, 2), B.block(0, 1, 2, 2).adj_);
}
