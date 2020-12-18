#include <test/unit/math/test_ad.hpp>
#include <cmath>
#include <complex>
#include <vector>
#include <type_traits>

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

namespace stan {
namespace test {
namespace internal {
template <typename T, require_not_eigen_t<T>* = nullptr>
auto test_abs(const T& x) {
  return stan::math::abs(x);
}
template <typename T, require_eigen_t<T>* = nullptr>
auto test_abs(const T& x) {
  return x.unaryExpr([](const auto x) { return stan::math::abs(x); });
}
}  // namespace internal
}  // namespace test
}  // namespace stan
TEST(mathMixMatFun, abs_varmat) {
  auto f = [](const auto& x1) {
    using stan::test::internal::test_abs;
    return test_abs(x1);
  };
  auto com_args = stan::test::internal::common_nonzero_args();
  std::vector<double> extra_args{0, -3, -2, 2, -17.3, -0.68, 2, 4};
  Eigen::VectorXd A(com_args.size() + extra_args.size());
  int i = 0;
  for (double x : com_args) {
    A(i) = x;
    ++i;
  }
  for (double x : extra_args) {
    A(i) = x;
    ++i;
  }
  stan::test::expect_ad_vector_matvar(f, A);
}
