#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, traceDot) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::trace_dot(x, y);
  };

  // 1x1
  Eigen::MatrixXd a11{{3}};
  Eigen::MatrixXd b11{{7}};
  stan::test::expect_ad(f, a11, b11);
  stan::test::expect_ad_matvar(f, a11, b11);

  // 0x0
  Eigen::MatrixXd m00(0, 0);
  stan::test::expect_ad(f, m00, m00);
  stan::test::expect_ad_matvar(f, m00, m00);

  // 2x2
  Eigen::MatrixXd a22{{1, 2}, {3, 4}};
  Eigen::MatrixXd b22{{5, 6}, {7, 8}};
  stan::test::expect_ad(f, a22, b22);
  stan::test::expect_ad_matvar(f, a22, b22);

  // 2x3 times 3x2 (rectangular)
  Eigen::MatrixXd a23{{1, 2, 3}, {4, 5, 6}};
  Eigen::MatrixXd b32{{7, 8}, {9, 10}, {11, 12}};
  stan::test::expect_ad(f, a23, b32);
  stan::test::expect_ad_matvar(f, a23, b32);

  // 3x2 times 2x3 (rectangular, transposed shape)
  stan::test::expect_ad(f, b32, a23);
  stan::test::expect_ad_matvar(f, b32, a23);

  // 3x3
  Eigen::MatrixXd a33{{1, -2, 3}, {0.5, 7, -1}, {2, 0, 4}};
  Eigen::MatrixXd b33{{3, 1, -2}, {0, 5, 1}, {-1, 2, 6}};
  stan::test::expect_ad(f, a33, b33);
  stan::test::expect_ad_matvar(f, a33, b33);

  // dimension mismatch: A.cols() != B.rows()
  stan::test::expect_ad(f, a22, b32);
  stan::test::expect_ad_matvar(f, a22, b32);

  stan::test::expect_ad(f, a23, b22);
  stan::test::expect_ad_matvar(f, a23, b22);

  // dimension mismatch: A.cols() == B.rows() but A.rows() != B.cols()
  stan::test::expect_ad(f, a23, a33);
  stan::test::expect_ad_matvar(f, a23, a33);
}
