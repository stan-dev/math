#include <test/unit/math/test_ad.hpp>
namespace rep_matrix_test { 
// y is scalar
auto f(int m, int n) {
  return [=](const auto& y) { return stan::math::rep_matrix(y, m, n); };
}

// y is row vector or column vector
auto g(int k) {
  return [=](const auto& y) { return stan::math::rep_matrix(y, k); };
}
}

TEST(MathMixMatFun, repMatrix) {
  double y = 3;
  stan::test::expect_ad(rep_matrix_test::f(0, 0), y);
  stan::test::expect_ad(rep_matrix_test::f(1, 1), y);
  stan::test::expect_ad(rep_matrix_test::f(2, 3), y);

  // illegal arguments---test throw
  stan::test::expect_ad(rep_matrix_test::f(-2, -1), y);

  Eigen::VectorXd a(3);
  a << 3, 3, 3;
  stan::test::expect_ad(rep_matrix_test::g(2), a);

  Eigen::RowVectorXd b(2);
  b << 2, 2;
  stan::test::expect_ad(rep_matrix_test::g(3), b);
}
