#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(MathMixMatFun, col) {
  auto f = [](int j) {
    return [=](const auto& y) { return stan::math::col(y, j); };
  };
  // normal and exception returns
  Eigen::MatrixXd a(0, 0);

  Eigen::MatrixXd b(1, 1);
  b << 1;

  Eigen::MatrixXd c(2, 2);
  c << 1, 2, 3, 4;

  Eigen::MatrixXd d(3, 2);
  d << 1, 2, 3, 4, 5, 6;

  Eigen::MatrixXd e(2, 3);
  e << 1, 2, 3, 4, 5, 6;

  for (const auto& x : std::vector<Eigen::MatrixXd>{a, b, c, d, e})
    for (size_t i = 0; i < 5; ++i)
      stan::test::expect_ad(f(i), x);
}
