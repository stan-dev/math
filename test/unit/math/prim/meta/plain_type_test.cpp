#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>

TEST(MetaTraitsPrimMat, plain_type_non_eigen) {
  test::expect_same_type<double, stan::plain_type_t<double>>();
  test::expect_same_type<std::vector<double>,
                         stan::plain_type_t<std::vector<double>>>();
  test::expect_same_type<stan::math::var,
                         stan::plain_type_t<stan::math::var>>();
  test::expect_same_type<std::vector<stan::math::var>,
                         stan::plain_type_t<std::vector<stan::math::var>>>();
}

TEST(MetaTraitsPrimMat, plain_type_eigen) {
  test::expect_same_type<Eigen::MatrixXd,
                         stan::plain_type_t<Eigen::MatrixXd>>();
  test::expect_same_type<Eigen::VectorXd,
                         stan::plain_type_t<Eigen::VectorXd>>();
  test::expect_same_type<Eigen::ArrayXXd,
                         stan::plain_type_t<Eigen::ArrayXXd>>();
  test::expect_same_type<stan::math::matrix_v,
                         stan::plain_type_t<stan::math::matrix_v>>();

  Eigen::MatrixXd m1, m2;
  Eigen::VectorXd v1, v2;
  Eigen::ArrayXXd a1, a2;
  stan::math::matrix_v mv1, mv2;
  test::expect_same_type<Eigen::MatrixXd,
                         stan::plain_type_t<decltype(m1 - m2)>>();
  test::expect_same_type<Eigen::VectorXd,
                         stan::plain_type_t<decltype(v1 + v2)>>();
  test::expect_same_type<Eigen::ArrayXXd,
                         stan::plain_type_t<decltype(a1 * a2.sin())>>();
  test::expect_same_type<stan::math::matrix_v,
                         stan::plain_type_t<decltype(mv1 * mv2)>>();
}
