#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mixScalFun, hypergeometric_2F1) {
  auto f = [](const auto& a, const auto& b, const auto& z) {
    using stan::math::hypergeometric_2F1;
    return hypergeometric_2F1(a[0], a[1], b[0], z);
  };

  stan::test::expect_ad(f, std::vector<double>{1, 1},
                           std::vector<double>{1},
                           0.6);
  stan::test::expect_ad(f, std::vector<double>{1, 31},
                           std::vector<double>{41},
                           0.5);
  stan::test::expect_ad(f, std::vector<double>{1, -2.1},
                           std::vector<double>{41},
                           0.8);
  stan::test::expect_ad(f, std::vector<double>{1, -1.6},
                           std::vector<double>{10.6},
                           -0.8);
  stan::test::expect_ad(f, std::vector<double>{-3.1, -2.2},
                           std::vector<double>{10.6},
                           0.3);
  stan::test::expect_ad(f, std::vector<double>{3.70975, 1},
                           std::vector<double>{2.70975},
                           -0.2);
  stan::test::expect_ad(f, std::vector<double>{3.70975, 1},
                           std::vector<double>{2.70975},
                           0.999696);
  stan::test::expect_ad(f, std::vector<double>{1, 12},
                           std::vector<double>{10},
                           1.0);
  stan::test::expect_ad(f, std::vector<double>{1, 12},
                           std::vector<double>{20},
                           1.2);
  stan::test::expect_ad(f, std::vector<double>{1, -0.6},
                           std::vector<double>{10.6},
                           0.3);
}

TEST(mixScalFun, hypergeometric_2F2) {
  auto f = [](const auto& a, const auto& b, const auto& z) {
    using stan::math::hypergeometric_pFq;
    return hypergeometric_pFq(a, b, z);
  };

  stan::test::expect_ad(f, std::vector<double>{4, 2},
                           std::vector<double>{6, 3},
                           4.0);
}
