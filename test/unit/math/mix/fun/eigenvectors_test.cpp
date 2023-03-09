#include <test/unit/math/test_ad.hpp>
#include <stdexcept>

TEST(mathMixFun, eigenvectors) {
  auto f = [](const auto& x) {
    using stan::math::eigenvectors;
    return eigenvectors(x);
  };
  for (const auto& x : stan::test::square_test_matrices(0, 2)) {
    stan::test::expect_ad(f, x);
  }

  Eigen::MatrixXd a32(3, 2);
  a32 << 3, -5, 7, -7.2, 9.1, -6.3;
  EXPECT_THROW(f(a32), std::invalid_argument);
}

TEST(mathMixFun, eigenvectorsComplex) {
  auto f = [](const auto& x) {
    using stan::math::eigenvectors;
    return eigenvectors(stan::math::to_complex(x, 0));
  };
  for (const auto& x : stan::test::square_test_matrices(0, 2)) {
    stan::test::expect_ad(f, x);
  }

  Eigen::MatrixXd a32(3, 2);
  a32 << 3, -5, 7, -7.2, 9.1, -6.3;
  EXPECT_THROW(f(a32), std::invalid_argument);
}
