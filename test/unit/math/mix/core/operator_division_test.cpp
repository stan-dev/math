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
