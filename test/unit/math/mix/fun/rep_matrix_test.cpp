#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, repMatrix) {
  // y is scalar
  auto f = [](int m, int n) {
    return [=](const auto& y) { return stan::math::rep_matrix(y, m, n); };
  };

  // y is row vector or column vector
  auto g = [](int k) {
    return [=](const auto& y) { return stan::math::rep_matrix(y, k); };
  };

  double y = 3;
  stan::test::expect_ad(f(0, 0), y);
  stan::test::expect_ad(f(1, 1), y);
  stan::test::expect_ad(f(2, 3), y);

  // illegal arguments---test throw
  stan::test::expect_ad(f(-2, -1), y);

  Eigen::VectorXd a(3);
  a << 3, 3, 3;
  stan::test::expect_ad(g(2), a);

  Eigen::RowVectorXd b(2);
  b << 2, 2;
  stan::test::expect_ad(g(3), b);
}

TEST(MathMixMatFun, repVarMatrix) {
  using stan::math::rep_matrix;
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  auto x_var = var(1.0);
  auto x
      = rep_matrix<var, var_value<Eigen::Matrix<double, -1, -1>>>(x_var, 5, 5);
  auto x_sum = sum(x);
  x_sum.grad();

  std::cout << "\nx_sum: \n"
            << "(val): \n"
            << x_sum.val() << "\n(adj): \n"
            << x_sum.adj() << "\n";
  std::cout << "\nx: \n"
            << "(val): \n"
            << x.val() << "\n(adj): \n"
            << x.adj() << "\n";
  std::cout << "\nx_var: \n"
            << "(val): \n"
            << x_var.val() << "\n(adj): \n"
            << x_var.adj() << "\n";
}
