#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/util.hpp>

TEST_F(mathMix, diagMatrix) {
  auto f = [](const auto& y) { return stan::math::diag_matrix(y); };

  for (int i = 0; i < 5; ++i) {
    Eigen::VectorXd a(i);
    for (int n = 0; n < a.size(); ++n)
      a(n) = n;
    stan::test::expect_ad(f, a);
  }
}
