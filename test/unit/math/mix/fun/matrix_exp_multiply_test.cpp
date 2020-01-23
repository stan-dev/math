#include <test/unit/math/test_ad.hpp>

auto f_idx(int i, int j) {
  return [=](const auto& x, const auto& y) {
    return stan::math::matrix_exp_multiply(x, y)(i, j);
  };
}

void expect_matrix_exp_multiply(int m, int n) {
  auto f = [](const auto& x, auto& y) {
    return stan::math::matrix_exp_multiply(x, y);
  };
  Eigen::MatrixXd a = Eigen::MatrixXd::Random(n, n);
  Eigen::MatrixXd b = Eigen::MatrixXd::Random(n, m);
  stan::test::expect_ad(f, a, b);
}

TEST(MathMixMatFun, matrixExpMultiply) {
  auto f = [](const auto& x, auto& y) {
    return stan::math::matrix_exp_multiply(x, y);
  };

  Eigen::MatrixXd m00(0, 0);
  stan::test::expect_ad(f, m00, m00);

  Eigen::MatrixXd m02(0, 2);
  stan::test::expect_ad(f, m00, m02);

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;

  Eigen::MatrixXd a(5, 5);
  a << -0.96871, 0.398827, 0.241306, 0.741373, 0.108926, 0.888077, -0.915624,
      -0.373344, 0.255238, 0.717304, -0.0899219, -0.898862, -0.800546,
      -0.222652, -0.271382, 0.683227, 0.827031, -0.780702, -0.104228, 0.885106,
      -0.996585, -0.097802, 0.739617, 0.235266, -0.0247717;
  Eigen::MatrixXd b(5, 5);
  b << -0.96871, 0.398827, 0.241306, 0.741373, 0.108926, 0.888077, -0.915624,
      -0.373344, 0.255238, 0.717304, -0.0899219, -0.898862, -0.800546,
      -0.222652, -0.271382, 0.683227, 0.827031, -0.780702, -0.104228, 0.885106,
      -0.996585, -0.097802, 0.739617, 0.235266, -0.0247717;
  // this matches previous tests in only testing the first position
  stan::test::expect_ad(tols, f_idx(0, 0), a, b);

  Eigen::MatrixXd c(1, 1);
  c << 1.3;
  Eigen::MatrixXd d(1, 1);
  d << -2.4;
  stan::test::expect_ad(f, c, d);

  Eigen::MatrixXd g(2, 2);
  Eigen::MatrixXd h(2, 2);
  g << 0.27, 0.73, 0.50, 0.50;
  h << 1, 2, 3, 4;
  stan::test::expect_ad(tols, f, g, h);

  Eigen::MatrixXd u(3, 3);
  Eigen::MatrixXd v(3, 3);
  u << -0.96, 0.39, 0.24, 0.74, 0.11, 0.88, -0.92, -0.37, 0.26;
  v << 0.13, 0.56, -0.46, 0.13, -0.11, -0.76, -0.23, -0.01, 0.37;
  stan::test::expect_ad(tols, f, u, v);

  // scaled down to size 3 from previous tests to test all args
  // at all orders
  Eigen::MatrixXd aa(3, 3);
  aa << -0.96, 0.39, 0.24, 0.74, 0.11, 0.88, -0.92, -0.37, 0.26;
  Eigen::MatrixXd bb(3, 1);
  bb << 1, 2, 3;
  stan::test::expect_ad(tols, f, aa, bb);

  Eigen::MatrixXd cc(1, 1);
  cc << 1.5;
  Eigen::MatrixXd dd(1, 3);
  dd << -3, -2, -1;
  stan::test::expect_ad(tols, f, cc, dd);
}
