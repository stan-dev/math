#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>

TEST(MetaTraitsPrimMat, plain_type_non_eigen) {
  EXPECT_SAME_TYPE(double, stan::plain_type_t<double>);
  EXPECT_SAME_TYPE(std::vector<double>,
                   stan::plain_type_t<std::vector<double>>);
  EXPECT_SAME_TYPE(stan::math::var, stan::plain_type_t<stan::math::var>);
  EXPECT_SAME_TYPE(std::vector<stan::math::var>,
                   stan::plain_type_t<std::vector<stan::math::var>>);
}

TEST(MetaTraitsPrimMat, plain_type_eigen) {
  EXPECT_SAME_TYPE(Eigen::MatrixXd, stan::plain_type_t<Eigen::MatrixXd>);
  EXPECT_SAME_TYPE(Eigen::VectorXd, stan::plain_type_t<Eigen::VectorXd>);
  EXPECT_SAME_TYPE(Eigen::ArrayXXd, stan::plain_type_t<Eigen::ArrayXXd>);
  EXPECT_SAME_TYPE(stan::math::matrix_v,
                   stan::plain_type_t<stan::math::matrix_v>);

  Eigen::MatrixXd m1, m2;
  Eigen::VectorXd v1, v2;
  Eigen::ArrayXXd a1, a2;
  stan::math::matrix_v mv1, mv2;
  EXPECT_SAME_TYPE(Eigen::MatrixXd, stan::plain_type_t<decltype(m1 - m2)>);
  EXPECT_SAME_TYPE(Eigen::VectorXd, stan::plain_type_t<decltype(v1 + v2)>);
  EXPECT_SAME_TYPE(Eigen::ArrayXXd,
                   stan::plain_type_t<decltype(a1 * a2.sin())>);
  EXPECT_SAME_TYPE(stan::math::matrix_v,
                   stan::plain_type_t<decltype(mv1 * mv2)>);
}
