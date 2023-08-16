#include <test/unit/math/test_ad.hpp>
#include <stdexcept>

TEST(mathMixFun, eigendecompose) {
  auto f = [](const auto& x) {
    using stan::math::eigendecompose;
    return std::get<0>(eigendecompose(x));
  };
  auto g = [](const auto& x) {
    using stan::math::eigendecompose;
    return std::get<1>(eigendecompose(x));
  };
  for (const auto& x : stan::test::square_test_matrices(0, 2)) {
    stan::test::expect_ad(f, x);
    stan::test::expect_ad(g, x);
  }

  Eigen::MatrixXd a32(3, 2);
  a32 << 3, -5, 7, -7.2, 9.1, -6.3;
  EXPECT_THROW(f(a32), std::invalid_argument);
  EXPECT_THROW(g(a32), std::invalid_argument);
}

TEST(mathMixFun, eigendecomposeComplex) {
  auto f = [](const auto& x) {
    using stan::math::eigendecompose;
    return std::get<0>(eigendecompose(stan::math::to_complex(x, 0)));
  };
  auto g = [](const auto& x) {
    using stan::math::eigendecompose;
    return std::get<1>(eigendecompose(stan::math::to_complex(x, 0)));
  };
  for (const auto& x : stan::test::square_test_matrices(0, 2)) {
    stan::test::expect_ad(f, x);
    stan::test::expect_ad(g, x);
  }

  Eigen::MatrixXd a32(3, 2);
  a32 << 3, -5, 7, -7.2, 9.1, -6.3;
  EXPECT_THROW(f(a32), std::invalid_argument);
  EXPECT_THROW(g(a32), std::invalid_argument);
}
