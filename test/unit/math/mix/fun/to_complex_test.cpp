#include <test/unit/math/test_ad.hpp>
#include <complex>

TEST(mathMixMatFun, to_complex) {
  auto f0 = []() { return stan::math::to_complex(); };
  auto f1 = [](const auto& x) { return stan::math::to_complex(x); };
  auto f2 = [](const auto& x, const auto& y) {
    return stan::math::to_complex(x, y);
  };

  EXPECT_EQ(std::complex<double>(), stan::math::to_complex());
  stan::test::expect_common_unary(f1);
  stan::test::expect_common_binary(f2);
}

TEST(mathMixMatFun, to_complex_vectorized) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::to_complex(x, y);
  };

  Eigen::VectorXd in1(2);
  in1 << 0.5, 3.4;
  std::vector<int> std_in2{3, 1};
  stan::test::expect_ad_vectorized_binary(f, in1, std_in2);

  Eigen::MatrixXd mat_in1 = in1.replicate(1, 2);
  std::vector<std::vector<int>> std_std_in2{std_in2, std_in2};
  stan::test::expect_ad_vectorized_binary(f, mat_in1, std_std_in2);
}
