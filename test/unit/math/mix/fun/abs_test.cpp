#include <test/unit/math/test_ad.hpp>
#include <cmath>
#include <complex>
#include <vector>
#include <type_traits>

TEST(mixFun, absBasics) {
  using stan::math::abs;
  int a = abs(1);

  double b = abs(-2.3);

  Eigen::Matrix<double, -1, -1> x(2, 3);
  x << 1, 2, 3, 4, 5, 6;
  Eigen::Matrix<double, -1, -1> y = abs(x);

  std::vector<int> u{1, 2, 3, 4};
  std::vector<int> v = abs(u);
}

TEST(mixFun, abs) {
  auto f = [](const auto& x) {
    using std::abs;
    return abs(x);
  };
  stan::test::expect_common_nonzero_unary(f);
  stan::test::expect_value(f, 0);
  stan::test::expect_value(f, 0.0);

  stan::test::expect_ad(f, -3);
  stan::test::expect_ad(f, -2);
  stan::test::expect_ad(f, 2);

  stan::test::expect_ad(f, -17.3);
  stan::test::expect_ad(f, -0.68);
  stan::test::expect_ad(f, 0.68);
  stan::test::expect_ad(f, 2.0);
  stan::test::expect_ad(f, 4.0);

  // not differentiable at zero
  for (double re : std::vector<double>{-4, -2.5, -1.5, -0.3, 1.3, 2.1, 3.9}) {
    for (double im : std::vector<double>{-4, -2.5, -1.5, -0.3, 1.3, 2.1, 3.9}) {
      stan::test::expect_ad(f, std::complex<double>(re, im));
    }
  }
}
TEST(mixFun, absReturnType) {
  // validate return types not overpromoted to complex by assignability
  std::complex<stan::math::var> a = 3;
  stan::math::var b = abs(a);

  std::complex<stan::math::fvar<double>> c = 3;
  stan::math::fvar<double> d = abs(c);
  SUCCEED();
}

TEST(mathMixMatFun, abs_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_nonzero_args;
  auto f = [](const auto& x1) {
    using stan::math::abs;
    return abs(x1);
  };
  std::vector<double> com_args = common_nonzero_args();
  std::vector<double> args{-3, 2, -0.68, 1};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}
