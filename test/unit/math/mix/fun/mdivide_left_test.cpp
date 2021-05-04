#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(MathMixMatFun, mdivideLeft) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::mdivide_left(x, y);
  };

  // size zero inputs
  Eigen::MatrixXd m00(0, 0);
  Eigen::MatrixXd m02(0, 2);
  Eigen::VectorXd v0(0);
  stan::test::expect_ad(f, m00, m00);
  stan::test::expect_ad(f, m00, m02);
  stan::test::expect_ad(f, m00, v0);
  stan::test::expect_ad_matvar(f, m00, m00);
  stan::test::expect_ad_matvar(f, m00, m02);
  stan::test::expect_ad_matvar(f, m00, v0);

  Eigen::MatrixXd aa(1, 1);
  aa << 1;
  Eigen::MatrixXd bb(1, 1);
  bb << 2;
  stan::test::expect_ad(f, aa, bb);
  stan::test::expect_ad_matvar(f, aa, bb);
  Eigen::MatrixXd b0(1, 0);
  stan::test::expect_ad(f, aa, b0);
  stan::test::expect_ad_matvar(f, aa, b0);

  Eigen::VectorXd cc(1);
  cc << 3;
  stan::test::expect_ad(f, aa, cc);
  stan::test::expect_ad_matvar(f, aa, cc);

  Eigen::MatrixXd a(2, 2);
  a << 2, 3, 3, 7;

  Eigen::MatrixXd b(2, 2);
  b << 3, 0, 0, 4;

  Eigen::MatrixXd c(2, 2);
  c << 12, 13, 15, 17;

  Eigen::MatrixXd d(2, 2);
  d << 2, 3, 5, 7;

  Eigen::MatrixXd e(2, 0);

  Eigen::VectorXd g(2);
  g << 12, 13;

  // matrix, matrix
  for (const auto& m1 : std::vector<Eigen::MatrixXd>{a, b, c, d}) {
    for (const auto& m2 : std::vector<Eigen::MatrixXd>{a, b, c, d, e}) {
      stan::test::expect_ad(f, m1, m2);
      stan::test::expect_ad_matvar(f, m1, m2);
    }
  }

  // matrix, vector
  for (const auto& m : std::vector<Eigen::MatrixXd>{a, b, c, d}) {
    stan::test::expect_ad(f, m, g);
    stan::test::expect_ad_matvar(f, m, g);
  }

  Eigen::MatrixXd v(5, 5);
  v << 20, 8, -9, 7, 5, 8, 20, 0, 4, 4, -9, 0, 20, 2, 5, 7, 4, 2, 20, -5, 5, 4,
      5, -5, 20;
  Eigen::VectorXd u(5);
  u << 62, 84, 84, 76, 108;
  stan::test::expect_ad(f, v, u);
  stan::test::expect_ad_matvar(f, v, u);

  Eigen::MatrixXd m33 = Eigen::MatrixXd::Zero(3, 3);
  Eigen::MatrixXd m44 = Eigen::MatrixXd::Zero(4, 4);
  Eigen::VectorXd v3 = Eigen::VectorXd::Zero(3);
  Eigen::VectorXd v4 = Eigen::VectorXd::Zero(4);
  Eigen::RowVectorXd rv3 = Eigen::RowVectorXd::Zero(3);
  Eigen::RowVectorXd rv4 = Eigen::RowVectorXd::Zero(4);

  // exceptions: wrong sizes
  stan::test::expect_ad(f, m33, m44);
  stan::test::expect_ad(f, m33, v4);
  stan::test::expect_ad_matvar(f, m33, m44);
  stan::test::expect_ad_matvar(f, m33, v4);

  // exceptions: wrong types
  stan::test::expect_ad(f, m33, rv3);
  stan::test::expect_ad_matvar(f, m33, rv3);
}
