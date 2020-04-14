#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

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

  stan::test::expect_ad(tols, f, p_rv, ds);
  stan::test::expect_ad(tols, f, p_rv, ds_v);
  stan::test::expect_ad(tols, f, p_rv, ds_rv);
}

TEST(mathMixMatFun, logMix) {
  std::vector<double> v1{-1.3};
  std::vector<double> v2{-1.3, 3.2};

  // old tests translated
  std::vector<double> a5{0.235, 0.152, 0.359, 0.254};
  std::vector<double> a7{0.514, 0.284, 0.112, 0.090};

  std::vector<double> i{-5.581, -6.254, -3.987, -10.221};
  std::vector<double> j{-6.785, -4.351, -5.847, -7.362};
  std::vector<double> k{-7.251, -10.510, -12.302, -3.587};
  std::vector<double> m{-1, -2, -3, -4};
  std::vector<double> n{-3.581, -8.114, -11.215, -5.658};

  expect_log_mix(a5, std::vector<std::vector<double>>{i, j, k});
  expect_log_mix(a7, std::vector<std::vector<double>>{i, m, n});

  // exception---mismatched sizes
  expect_log_mix(a5, std::vector<std::vector<double>>{v1, v2, j});
}
