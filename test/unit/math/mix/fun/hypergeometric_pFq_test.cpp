#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mixScalFun, hypergeometric_pFq_ad) {
  auto f = [](const auto& a, const auto& b, const auto& z) {
    using stan::math::hypergeometric_pFq;
    return hypergeometric_pFq(a, b, z);
  };

  stan::test::expect_ad(f, std::vector<double>{4, 2},
                           std::vector<double>{6, 3},
                           4.0);

  stan::test::expect_ad(f, std::vector<double>{2, 3},
                           std::vector<double>{2, 4, 5},
                           1.0);

  stan::test::expect_ad(f, std::vector<double>{1, 2, 3, 4},
                           std::vector<double>{5, 6, 7},
                           1.5);

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

  stan::test::expect_ad(f, std::vector<double>{1, 31, -27},
                           std::vector<double>{19, -41},
                           1.0);
}
