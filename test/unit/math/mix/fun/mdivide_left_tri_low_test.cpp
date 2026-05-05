#include <test/unit/math/test_ad.hpp>
#include <limits>

namespace {
auto f = [](const auto& x, const auto& y) {
  // no need to symmetrize because only uses view
  return stan::math::mdivide_left_tri_low(x, y);
};

TEST(MathMixMatFun, mdivideLeftTriLow1) {
  // size zero inputs
  Eigen::MatrixXd m00(0, 0);
  stan::test::expect_ad(f, m00, m00);
}
TEST(MathMixMatFun, mdivideLeftTriLow2) {
  // size zero inputs
  Eigen::MatrixXd m00(0, 0);
  Eigen::MatrixXd m02(0, 2);
  stan::test::expect_ad(f, m00, m02);
}
TEST(MathMixMatFun, mdivideLeftTriLow3) {
  // size zero inputs
  Eigen::MatrixXd m00(0, 0);
  Eigen::VectorXd v0(0);
  stan::test::expect_ad(f, m00, v0);
}
TEST(MathMixMatFun, mdivideLeftTriLow4) {
  // signature 1 of 2: matrix-matrix
  Eigen::MatrixXd aa(1, 1);
  aa << 1;
  Eigen::MatrixXd bb(1, 1);
  bb << 2;
  stan::test::expect_ad(f, aa, bb);
}
TEST(MathMixMatFun, mdivideLeftTriLow5) {
  // signature 2 of 2: matrix-vector
  Eigen::MatrixXd aa(1, 1);
  aa << 1;
  Eigen::VectorXd cc(1);
  cc << 3;
  stan::test::expect_ad(f, aa, cc);
}

TEST(MathMixMatFun, mdivideLeftTriLow6a) {
  Eigen::MatrixXd a(2, 2);
  a << 1, std::numeric_limits<double>::quiet_NaN(), -3, 5;
  Eigen::MatrixXd b(2, 2);
  b << 2, 5, 12, 109;
  auto ff = [&a](auto&& b) { return f(a, b); };
  stan::test::expect_ad(ff, b);
}

TEST(MathMixMatFun, mdivideLeftTriLow6b) {
  Eigen::MatrixXd a(2, 2);
  a << 1, std::numeric_limits<double>::quiet_NaN(), -3, 5;
  Eigen::MatrixXd b(2, 2);
  b << 2, 5, 12, 109;
  auto ff = [&b](auto&& a) { return f(a, b); };
  stan::test::expect_ad(ff, a);
}

TEST(MathMixMatFun, mdivideLeftTriLow6c) {
  Eigen::MatrixXd a(2, 2);
  a << 1, std::numeric_limits<double>::quiet_NaN(), -3, 5;
  Eigen::MatrixXd b(2, 2);
  b << 2, 5, 12, 109;
  stan::test::expect_ad(f, a, b);
}

TEST(MathMixMatFun, mdivideLeftTriLow11) {
  Eigen::MatrixXd a(2, 2);
  a << 1, std::numeric_limits<double>::quiet_NaN(), -3, 5;
  Eigen::VectorXd c(2);
  c << 3, 5;
  stan::test::expect_ad(f, a, c);
}

TEST(MathMixMatFun, mdivideLeftTriLow12a2) {
  Eigen::MatrixXd a2(2, 2);
  a2 << 1, 0, -3, 5;
  Eigen::MatrixXd b2(2, 2);
  b2 << 2, 5, 12, 109;
  auto ff = [&b2](auto&& a2) { return f(a2, b2); };
  stan::test::expect_ad(ff, a2);
}
TEST(MathMixMatFun, mdivideLeftTriLow12b2) {
  Eigen::MatrixXd a2(2, 2);
  a2 << 1, 0, -3, 5;
  Eigen::MatrixXd b2(2, 2);
  b2 << 2, 5, 12, 109;
  auto ff = [&a2](auto&& b2) { return f(a2, b2); };
  stan::test::expect_ad(ff, b2);
}
TEST(MathMixMatFun, mdivideLeftTriLow12) {
  Eigen::MatrixXd a2(2, 2);
  a2 << 1, 0, -3, 5;
  Eigen::MatrixXd b2(2, 2);
  b2 << 2, 5, 12, 109;
  stan::test::expect_ad(f, a2, b2);
}
TEST(MathMixMatFun, mdivideLeftTriLow13) {
  Eigen::MatrixXd y(3, 3);
  y << 1, 0, 0, 2, 3, 0, 4, 5, 6;
  Eigen::VectorXd z(3);
  z << 1, 2, 3;
  stan::test::expect_ad(f, y, z);
}
TEST(MathMixMatFun, mdivideLeftTriLow14) {
  Eigen::MatrixXd u(3, 3);
  u << 1, 2, 3, 6, 5, 4, 7, 8, 9;
  Eigen::MatrixXd y(3, 3);
  y << 1, 0, 0, 2, 3, 0, 4, 5, 6;
  stan::test::expect_ad(f, u, y);
}

TEST(MathMixMatFun, mdivideLeftTriLow15) {
  Eigen::MatrixXd b(2, 2);
  b << 2, 5, 12, 109;
  Eigen::MatrixXd d(3, 2);
  d << 1, 2, 5, -19, 73, 31;
  stan::test::expect_ad(f, d, b);
}
TEST(MathMixMatFun, mdivideLeftTriLow21) {
  Eigen::MatrixXd d(3, 2);
  d << 1, 2, 5, -19, 73, 31;
  Eigen::VectorXd c(2);
  c << 3, 5;
  stan::test::expect_ad(f, d, c);
}
TEST(MathMixMatFun, mdivideLeftTriLow22) {
  Eigen::MatrixXd a(2, 2);
  a << 1, std::numeric_limits<double>::quiet_NaN(), -3, 5;
  Eigen::VectorXd e(3);
  e << 1, 2, 9;
  stan::test::expect_ad(f, a, e);
}
TEST(MathMixMatFun, mdivideLeftTriLow23) {
  // exceptions: wrong sizes
  Eigen::MatrixXd m33 = Eigen::MatrixXd::Zero(3, 3);
  Eigen::MatrixXd m44 = Eigen::MatrixXd::Zero(4, 4);
  stan::test::expect_ad(f, m33, m44);
}
TEST(MathMixMatFun, mdivideLeftTriLow24) {
  Eigen::MatrixXd m33 = Eigen::MatrixXd::Zero(3, 3);
  Eigen::VectorXd v4 = Eigen::VectorXd::Zero(4);
  stan::test::expect_ad(f, m33, v4);
}
TEST(MathMixMatFun, mdivideLeftTriLow25) {
  // exceptions: wrong types
  Eigen::MatrixXd m33 = Eigen::MatrixXd::Zero(3, 3);
  Eigen::RowVectorXd rv3 = Eigen::RowVectorXd::Zero(3);
  stan::test::expect_ad(f, m33, rv3);
}

}  // namespace
