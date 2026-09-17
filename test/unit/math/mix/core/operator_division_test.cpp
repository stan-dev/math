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
