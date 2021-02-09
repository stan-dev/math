#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

TEST(mathMixMatFun, dotProduct) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::dot_product(x, y);
  };

  // size 0
  Eigen::VectorXd v0;
  Eigen::RowVectorXd rv0;
  std::vector<double> sv0;
  stan::test::expect_ad(f, v0, v0);
  stan::test::expect_ad(f, v0, rv0);
  stan::test::expect_ad(f, rv0, v0);
  stan::test::expect_ad(f, rv0, rv0);
  stan::test::expect_ad(f, sv0, sv0);

  stan::test::expect_ad_matvar(f, v0, v0);
  stan::test::expect_ad_matvar(f, v0, rv0);
  stan::test::expect_ad_matvar(f, rv0, v0);
  stan::test::expect_ad_matvar(f, rv0, rv0);

  // size 1
  Eigen::VectorXd v1(1);
  v1 << 1;
  Eigen::RowVectorXd rv1(1);
  rv1 << 1;
  std::vector<double> sv1{1};
  Eigen::VectorXd v1b(1);
  v1b << 2;
  Eigen::RowVectorXd rv1b(1);
  rv1b << 2;
  std::vector<double> sv1b{2};
  stan::test::expect_ad(f, v1, v1b);
  stan::test::expect_ad(f, rv1, v1b);
  stan::test::expect_ad(f, v1, rv1b);
  stan::test::expect_ad(f, rv1, rv1b);
  stan::test::expect_ad(f, sv1, sv1b);

  stan::test::expect_ad_matvar(f, v1, v1b);
  stan::test::expect_ad_matvar(f, rv1, v1b);
  stan::test::expect_ad_matvar(f, v1, rv1b);
  stan::test::expect_ad_matvar(f, rv1, rv1b);

  // size 2
  Eigen::VectorXd v2(2);
  v2 << 1, 2;
  Eigen::RowVectorXd rv2(2);
  rv2 << 1, 2;
  std::vector<double> sv2 = {1.0, 2.0};
  Eigen::VectorXd v2b(2);
  v2b << 10, 100;
  Eigen::RowVectorXd rv2b(2);
  rv2b << 10, 100;
  std::vector<double> sv2b = {10.0, 100.0};
  stan::test::expect_ad(f, v2, v2b);
  stan::test::expect_ad(f, rv2, v2b);
  stan::test::expect_ad(f, v2, rv2b);
  stan::test::expect_ad(f, rv2, rv2b);
  stan::test::expect_ad(f, sv2, sv2b);

  stan::test::expect_ad_matvar(f, v2, v2b);
  stan::test::expect_ad_matvar(f, rv2, v2b);
  stan::test::expect_ad_matvar(f, v2, rv2b);
  stan::test::expect_ad_matvar(f, rv2, rv2b);

  // size 3
  Eigen::VectorXd v3(3);
  v3 << 1, 3, -5;
  Eigen::RowVectorXd rv3(3);
  rv3 << 1, 3, -5;
  std::vector<double> sv3 = {1.0, 2.0, 3.0};
  Eigen::VectorXd v3b(3);
  v3b << 4, -2, -1;
  Eigen::RowVectorXd rv3b(3);
  rv3b << 4, -2, -1;
  std::vector<double> sv3b = {4.0, -2.0, -1.0};
  stan::test::expect_ad(f, v3, v3b);
  stan::test::expect_ad(f, v3, rv3b);
  stan::test::expect_ad(f, rv3, v3b);
  stan::test::expect_ad(f, rv3, rv3b);
  stan::test::expect_ad(f, sv3, sv3b);

  stan::test::expect_ad_matvar(f, v3, v3b);
  stan::test::expect_ad_matvar(f, v3, rv3b);
  stan::test::expect_ad_matvar(f, rv3, v3b);
  stan::test::expect_ad_matvar(f, rv3, rv3b);

  // size 3, another case (originally from rev)
  v3 << -1, 0, 1;
  rv3 << -1, 0, 1;
  sv3 = {-1.0, 0.0, 1.0};
  v3b << 1, 2, 3;
  rv3b << 1, 2, 3;
  sv3b = {1.0, 2.0, 3.0};
  stan::test::expect_ad(f, v3, v3b);
  stan::test::expect_ad(f, v3, rv3b);
  stan::test::expect_ad(f, rv3, v3b);
  stan::test::expect_ad(f, rv3, rv3b);
  stan::test::expect_ad(f, sv3, sv3b);

  stan::test::expect_ad_matvar(f, v3, v3b);
  stan::test::expect_ad_matvar(f, v3, rv3b);
  stan::test::expect_ad_matvar(f, rv3, v3b);
  stan::test::expect_ad_matvar(f, rv3, rv3b);

  // following throw due to mismatched size
  stan::test::expect_ad(f, v1, v3);
  stan::test::expect_ad(f, v1, rv3);
  stan::test::expect_ad(f, rv1, v3);
  stan::test::expect_ad(f, rv1, rv3);
  stan::test::expect_ad(f, sv1, sv3);

  stan::test::expect_ad_matvar(f, v1, v3);
  stan::test::expect_ad_matvar(f, v1, rv3);
  stan::test::expect_ad_matvar(f, rv1, v3);
  stan::test::expect_ad_matvar(f, rv1, rv3);

  // eliminated tests that matrix args throw in ad because
  // they don't compile in prim;  they shouldn't compile in ad, either
}
