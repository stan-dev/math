#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(MathMixMatFun, norm1) {
  auto f = [](const auto& y) { return stan::math::norm1(y); };

  Eigen::VectorXd x0(0);

  Eigen::VectorXd x1(1);
  x1 << 2;

  Eigen::VectorXd x2(2);
  x2 << 2, 3;

  Eigen::VectorXd x3(3);
  x3 << 2, 3, 4;

  for (const auto& a : std::vector<Eigen::VectorXd>{x0, x1, x2, x3}) {
    stan::test::expect_ad(f, a);
    stan::test::expect_ad_matvar(f, a);
  }

  std::vector<Eigen::VectorXd> vec_eig{x3, x3, x3};
  stan::test::expect_ad(f, vec_eig);
  stan::test::expect_ad_matvar(f, vec_eig);

  std::vector<double> vec_real{1, 2, 3, 4};
  stan::test::expect_ad(f, vec_real);
}
