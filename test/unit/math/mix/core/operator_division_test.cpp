#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(mathMixCore, operatorDivision) {
  auto f = [](const auto& x1, const auto& x2) { return x1 / x2; };
  bool disable_lhs_int = true;
  stan::test::expect_common_binary(f, disable_lhs_int);

  std::vector<double> common_finite = {-2.9, -1, -0.0, 0.0, 1, 1.39};
  std::vector<double> common_finite_nz = {-3.1, -1, 1, 2.7};
  for (auto re1 : common_finite) {
    for (auto im1 : common_finite_nz) {
      for (auto re2 : common_finite) {
        for (auto im2 : common_finite_nz) {
          stan::test::expect_ad(f, std::complex<double>(re1, im1),
                                std::complex<double>(re2, im2));
        }
      }
    }
  }
  for (auto re1 : common_finite) {
    for (auto re2 : common_finite) {
      for (auto im2 : common_finite_nz) {
        stan::test::expect_ad(f, re1, std::complex<double>(re2, im2));
      }
    }
  }
  for (auto re1 : common_finite) {
    for (auto im1 : common_finite_nz) {
      for (auto re2 : common_finite) {
        stan::test::expect_ad(f, std::complex<double>(re1, im1), re2);
      }
    }
  }
}

namespace stan {
namespace test {
struct operator_divide_tester {
  template <typename T1, typename T2,
            require_any_var_matrix_t<T1, T2>* = nullptr>
  auto operator()(const T1& x, const T2& y) const {
    return x / y;
  }
  template <typename T1, typename T2, require_any_eigen_t<T1, T2>* = nullptr,
            require_all_not_var_matrix_t<T1, T2>* = nullptr>
  auto operator()(const T1& x, const T2& y) const {
    return (stan::math::as_array_or_scalar(x)
            / stan::math::as_array_or_scalar(y))
        .matrix()
        .eval();
  }
};
}  // namespace test
}  // namespace stan

TEST(mathMixCore, operatorDivisionVarMat) {
  Eigen::MatrixXd mat1(2, 2);
  mat1 << -2, -1, 0.5, 2.8;
  Eigen::MatrixXd mat2 = mat1.reverse();
  stan::test::expect_ad_matvar(stan::test::operator_divide_tester{}, mat1,
                               mat2);
  stan::test::expect_ad_matvar(stan::test::operator_divide_tester{}, mat1, 2.0);
  stan::test::expect_ad_matvar(stan::test::operator_divide_tester{}, 2.0, mat2);
}

TEST(mathMixCore, division_extreme_tangents) {
  using namespace stan::math;
  const auto check
      = [](double a, double b, double da, double db, double expected) {
          fvar<double> x(a, da), y(b, db);
          const auto result = x / y;
          EXPECT_NEAR(expected, result.d_, std::abs(expected) * 1e-12 + 1e-323);
          x /= y;
          EXPECT_EQ(result.val_, x.val_);
          EXPECT_EQ(result.d_, x.d_);
        };
  check(1, 1e160, 0, 1, -1e-320);
  check(1e308, 1e308, 2, 2, 0);
  // The quotient underflows, but its product with the tangent is finite.
  check(1e-200, 1e150, 0, 1e308, -1e-192);
  check(1e-200, 1e-150, 0, 1e-200, -1e-100);
  check(1e-300, 1e-100, 0, 1e-150, -1e-250);
  check(1, std::ldexp(1.0, -1024), 1, std::ldexp(1.0, -1024), 0);
  check(1, 1e-310, 0, 1e-320, -(1e-320 / 1e-310) / 1e-310);
  {
    nested_rev_autodiff nested;
    var x = 1e-200, y = 1e150;
    var result = 1e308 * (x / y);
    result.grad();
    EXPECT_NEAR(-1e-192, y.adj(), 1e-204);
  }
}
