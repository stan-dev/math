#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

void expect_log_mix(const std::vector<double>& p,
                    const std::vector<double>& d) {
  auto f
      = [](const auto& x, const auto& y) { return stan::math::log_mix(x, y); };
  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;

  auto p_v = stan::test::to_vector(p);
  auto p_rv = stan::test::to_row_vector(p);
  auto d_v = stan::test::to_vector(d);
  auto d_rv = stan::test::to_row_vector(d);

  stan::test::expect_ad(tols, f, p, d);
  stan::test::expect_ad(tols, f, p, d_v);
  stan::test::expect_ad(tols, f, p, d_rv);
  stan::test::expect_ad(tols, f, p_v, d);
  stan::test::expect_ad(tols, f, p_v, d_v);
  stan::test::expect_ad(tols, f, p_v, d_rv);
  stan::test::expect_ad(tols, f, p_rv, d);
  stan::test::expect_ad(tols, f, p_rv, d_v);
  stan::test::expect_ad(tols, f, p_rv, d_rv);
}

std::vector<Eigen::VectorXd> to_vectors(
    const std::vector<std::vector<double>>& xs) {
  std::vector<Eigen::VectorXd> ys;
  for (const auto& x : xs)
    ys.push_back(stan::test::to_vector(x));
  return ys;
}
std::vector<Eigen::RowVectorXd> to_row_vectors(
    const std::vector<std::vector<double>>& xs) {
  std::vector<Eigen::RowVectorXd> ys;
  for (const auto& x : xs)
    ys.push_back(stan::test::to_row_vector(x));
  return ys;
}

void expect_log_mix(const std::vector<double>& p,
                    const std::vector<std::vector<double>>& ds) {
  auto f
      = [](const auto& x, const auto& y) { return stan::math::log_mix(x, y); };
  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;

  auto p_v = stan::test::to_vector(p);
  auto p_rv = stan::test::to_row_vector(p);

  auto ds_v = to_vectors(ds);
  auto ds_rv = to_row_vectors(ds);

  stan::test::expect_ad(tols, f, p, ds);
  stan::test::expect_ad(tols, f, p, ds_v);
  stan::test::expect_ad(tols, f, p, ds_rv);
  stan::test::expect_ad(tols, f, p_v, ds);
  stan::test::expect_ad(tols, f, p_v, ds_v);
  stan::test::expect_ad(tols, f, p_v, ds_rv);
  stan::test::expect_ad(tols, f, p_rv, ds);
  stan::test::expect_ad(tols, f, p_rv, ds_v);
  stan::test::expect_ad(tols, f, p_rv, ds_rv);
}

TEST(mathMixMatFun, logMix) {
  std::vector<double> u1{0.999};  // using 1 overflows with finite diffs
  std::vector<double> v1{-1.3};
  expect_log_mix(u1, v1);

  std::vector<double> u2{0.3, 0.7};
  std::vector<double> v2{-1.3, 3.2};
  expect_log_mix(u2, v2);

  // size exceptions
  expect_log_mix(u1, v2);
  expect_log_mix(u2, v1);

  // old tests translated
  std::vector<double> a{0.112, 0.214, 0.305, 0.369};
  std::vector<double> a2{0.15, 0.20, 0.40, 0.25};
  std::vector<double> a3{0.13, 0.22, 0.38, 0.27};
  std::vector<double> a4{0.03, 0.21, 0.63, 0.13};
  std::vector<double> a5{0.235, 0.152, 0.359, 0.254};
  std::vector<double> a6{0.15, 0.70, 0.10, 0.05};
  std::vector<double> a7{0.514, 0.284, 0.112, 0.090};

  std::vector<double> b{-5.983, -11.215, -6.836, -5.438};
  std::vector<double> c{-10.365, -12.443, -15.091, -19.115};
  std::vector<double> d{-4.117, -8.132, -7.931, -12.115};
  std::vector<double> e{-2.15, -3.89, -2.18, -8.82};
  std::vector<double> g{-3.15, -0.21, -10.55, -7.24};
  std::vector<double> h{-19.41, -8.14, -2.18, -9.13};
  std::vector<double> i{-5.581, -6.254, -3.987, -10.221};
  std::vector<double> j{-6.785, -4.351, -5.847, -7.362};
  std::vector<double> k{-7.251, -10.510, -12.302, -3.587};
  std::vector<double> m{-1, -2, -3, -4};
  std::vector<double> n{-3.581, -8.114, -11.215, -5.658};

  expect_log_mix(a, b);
  expect_log_mix(a, c);
  expect_log_mix(a, d);
  expect_log_mix(a2, e);
  expect_log_mix(a3, g);
  expect_log_mix(a4, h);
  expect_log_mix(a5, i);
  expect_log_mix(a6, m);
  expect_log_mix(a7, n);

  expect_log_mix(a5, std::vector<std::vector<double>>{i, j, k});
  expect_log_mix(a7, std::vector<std::vector<double>>{i, m, n});

  // exception---mismatched sizes
  expect_log_mix(a5, std::vector<std::vector<double>>{v1, v2, j});
}
