#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, mdivideRightSpd) {
  auto f = [](const auto& x, const auto& y) {
    if (y.rows() != y.cols())
      return stan::math::mdivide_right_spd(x, y);
    auto y_sym = ((y + y.transpose()) * 0.5).eval();  // sym for finite diffs
    return stan::math::mdivide_right_spd(x, y_sym);
  };

  // size zero inputs
  Eigen::MatrixXd m00(0, 0);
  Eigen::MatrixXd m20(2, 0);
  Eigen::RowVectorXd rv0(0);
  stan::test::expect_ad(f, m00, m00);
  stan::test::expect_ad(f, m20, m00);
  stan::test::expect_ad(f, rv0, m00);

  Eigen::MatrixXd aa(1, 1);
  aa << 1;
  Eigen::MatrixXd bb(1, 1);
  bb << 2;
  stan::test::expect_ad(f, aa, bb);
  Eigen::MatrixXd a0(0, 1);
  stan::test::expect_ad(f, a0, bb);

  Eigen::RowVectorXd cc(1);
  cc << 3;
  stan::test::expect_ad(f, cc, aa);

  Eigen::MatrixXd a(2, 2);
  a << 2, 3, 3, 7;

  Eigen::MatrixXd b(2, 2);
  b << 3, 0, 0, 4;

  Eigen::MatrixXd c(2, 2);
  c << 12, 13, 15, 17;

  Eigen::RowVectorXd d(2);
  d << 12, 13;

  // matrix / matrix : these have spd second args as required
  stan::test::expect_ad(f, a, a);
  stan::test::expect_ad(f, b, b);
  stan::test::expect_ad(f, c, a);
  stan::test::expect_ad(f, c, b);
  // row_vector / matrix: ditto
  stan::test::expect_ad(f, d, a);
  stan::test::expect_ad(f, d, b);

  Eigen::MatrixXd m33 = Eigen::MatrixXd::Zero(3, 3);
  Eigen::MatrixXd m44 = Eigen::MatrixXd::Zero(4, 4);
  Eigen::VectorXd v3 = Eigen::VectorXd::Zero(3);
  Eigen::RowVectorXd rv3 = Eigen::RowVectorXd::Zero(3);
  Eigen::RowVectorXd rv4 = Eigen::RowVectorXd::Zero(4);

  // exceptions: not symmetric
  stan::test::expect_ad(f, a, c);
  stan::test::expect_ad(f, d, c);

  // exceptions: not pos def
  stan::test::expect_ad(f, m33, m33);
  stan::test::expect_ad(f, rv3, m33);

  // exceptions: wrong sizes
  stan::test::expect_ad(f, m33, m44);
  stan::test::expect_ad(f, rv4, m33);

  // exceptions: wrong types
  stan::test::expect_ad(f, v3, m33);

  stan::math::recover_memory();
}
