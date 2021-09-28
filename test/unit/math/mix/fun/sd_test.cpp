#include <test/unit/math/test_ad.hpp>

template <typename T>
void expect_sd(const T& m, const stan::test::ad_tolerances& tols
                           = stan::test::ad_tolerances()) {
  auto f = [](const auto& x) { return stan::math::sd(x); };
  Eigen::VectorXd v(m.size());
  Eigen::RowVectorXd rv(m.size());
  std::vector<double> st(m.size());

  for (int i = 0; i < m.size(); ++i) {
    v(i) = m(i);
    rv(i) = m(i);
    st[i] = m(i);
  }

  std::vector<std::vector<double>> stst = {st, st};
  std::vector<Eigen::VectorXd> stv = {v, v};
  std::vector<Eigen::RowVectorXd> strv = {rv, rv};
  std::vector<Eigen::MatrixXd> stm = {m, m};

  stan::test::expect_ad(tols, f, v);
  stan::test::expect_ad(tols, f, rv);
  stan::test::expect_ad(tols, f, m);
  stan::test::expect_ad(tols, f, st);
  stan::test::expect_ad(tols, f, stv);
  stan::test::expect_ad(tols, f, strv);
  stan::test::expect_ad(tols, f, stm);
  stan::test::expect_ad(tols, f, stst);

  stan::test::expect_ad_matvar(tols, f, v);
  stan::test::expect_ad_matvar(tols, f, rv);
  stan::test::expect_ad_matvar(tols, f, m);
  stan::test::expect_ad_matvar(tols, f, stv);
  stan::test::expect_ad_matvar(tols, f, strv);
  stan::test::expect_ad_matvar(tols, f, stm);
}

TEST(MathMixMatFun, sd) {
  Eigen::MatrixXd a(0, 0);
  expect_sd(a);

  Eigen::MatrixXd b(1, 1);
  b << 1;
  expect_sd(b);

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-1;
  tols.hessian_fvar_hessian_ = 1e-1;
  Eigen::MatrixXd c(3, 1);
  c << 100, 0, -3;
  expect_sd(c, tols);

  Eigen::MatrixXd d(1, 3);
  d << 100, 0, -3;
  expect_sd(d, tols);

  Eigen::MatrixXd e(3, 2);
  e << -100, 0, 1, 20, -40, 2;
  expect_sd(e, tols);

  Eigen::MatrixXd g(6, 1);
  g << 1, 2, 3, 4, 5, 6;
  expect_sd(g, tols);
}
