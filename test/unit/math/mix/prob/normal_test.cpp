#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>

TEST_F(AgradRev, mathMixScalFun_normal_lpdf) {
  auto f = [](const double mu, const double sigma) {
    return [=](const auto& y) { return stan::math::normal_lpdf(y, mu, sigma); };
  };

  stan::test::expect_ad(f(0, 1), -2.3);
  stan::test::expect_ad(f(0, 1), 0.0);
  stan::test::expect_ad(f(0, 1), 1.7);
}

TEST_F(AgradRev, mathMixScalFun_normal_lpdf_vectors) {
  auto f = [](const auto& y, const auto& mu, const auto& sigma) {
    return stan::math::normal_lpdf(y, mu, sigma);
  };
  Eigen::VectorXd y(3);
  y << -0.3, 1.2, 2.5;
  Eigen::VectorXd mu(3);
  mu << 0.1, -0.4, 2.0;
  Eigen::VectorXd sigma(3);
  sigma << 0.5, 1.3, 2.2;
  stan::test::expect_ad(f, y, mu, sigma);
  stan::test::expect_ad(f, y, mu, 1.7);
  stan::test::expect_ad(f, y, 0.2, sigma);
  stan::test::expect_ad(f, 0.4, mu, sigma);
  stan::test::expect_ad(f, y, 0.2, 1.7);
  stan::test::expect_ad(f, 0.4, 0.2, sigma);
}

namespace normal_lpdf_test {
inline void append_adjoints(std::vector<double>& out,
                            const stan::math::var& x) {
  out.push_back(x.adj());
}

template <typename T, stan::require_var_matrix_t<T>* = nullptr>
void append_adjoints(std::vector<double>& out, const T& x) {
  for (Eigen::Index i = 0; i < x.size(); ++i) {
    out.push_back(x.adj().coeff(i));
  }
}

template <typename T, stan::require_eigen_vt<stan::is_var, T>* = nullptr>
void append_adjoints(std::vector<double>& out, const T& x) {
  for (Eigen::Index i = 0; i < x.size(); ++i) {
    out.push_back(x.coeff(i).adj());
  }
}

template <typename T, stan::require_st_arithmetic<T>* = nullptr>
void append_adjoints(std::vector<double>& out, const T& /* x */) {}

// Value and adjoints of normal_lpdf<propto>(y, mu, sigma) with the
// arguments made by make_y, make_mu and make_sigma.
template <bool propto, typename MakeY, typename MakeMu, typename MakeSigma>
std::vector<double> value_and_adjoints(const MakeY& make_y,
                                       const MakeMu& make_mu,
                                       const MakeSigma& make_sigma) {
  auto y = make_y();
  auto mu = make_mu();
  auto sigma = make_sigma();
  stan::math::var lp = stan::math::normal_lpdf<propto>(y, mu, sigma);
  lp.grad();
  std::vector<double> out{lp.val()};
  append_adjoints(out, y);
  append_adjoints(out, mu);
  append_adjoints(out, sigma);
  stan::math::recover_memory();
  return out;
}

template <bool propto>
void expect_varmat_matches_matrix_var() {
  using stan::math::var;
  using stan::math::var_value;
  using vm = var_value<Eigen::VectorXd>;
  Eigen::VectorXd y(4);
  y << -0.3, 1.2, 2.5, 0.0;
  Eigen::VectorXd mu(4);
  mu << 0.1, -0.4, 2.0, 0.7;
  Eigen::VectorXd sigma(4);
  sigma << 0.5, 1.3, 2.2, 0.9;
  auto data_y = [&]() { return y; };
  auto scalar_sigma = [&]() { return var(1.7); };
  auto expect_same
      = [](const std::vector<double>& a, const std::vector<double>& b) {
          ASSERT_EQ(a.size(), b.size());
          for (size_t i = 0; i < a.size(); ++i) {
            EXPECT_FLOAT_EQ(a[i], b[i]) << "entry " << i;
          }
        };
  // mu autodiff vector, sigma scalar
  expect_same(
      value_and_adjoints<propto>(
          data_y, [&]() { return stan::math::to_var(mu); }, scalar_sigma),
      value_and_adjoints<propto>(
          data_y, [&]() { return vm(mu); }, scalar_sigma));
  // mu and sigma autodiff vectors
  expect_same(
      value_and_adjoints<propto>(
          data_y, [&]() { return stan::math::to_var(mu); },
          [&]() { return stan::math::to_var(sigma); }),
      value_and_adjoints<propto>(
          data_y, [&]() { return vm(mu); }, [&]() { return vm(sigma); }));
  // y autodiff vector, mu data vector, sigma scalar
  expect_same(
      value_and_adjoints<propto>([&]() { return stan::math::to_var(y); },
                                 [&]() { return mu; }, scalar_sigma),
      value_and_adjoints<propto>([&]() { return vm(y); }, [&]() { return mu; },
                                 scalar_sigma));
  // all autodiff vectors
  expect_same(
      value_and_adjoints<propto>([&]() { return stan::math::to_var(y); },
                                 [&]() { return stan::math::to_var(mu); },
                                 [&]() { return stan::math::to_var(sigma); }),
      value_and_adjoints<propto>([&]() { return vm(y); },
                                 [&]() { return vm(mu); },
                                 [&]() { return vm(sigma); }));
}
}  // namespace normal_lpdf_test

TEST_F(AgradRev, mathMixScalFun_normal_lpdf_varmat) {
  normal_lpdf_test::expect_varmat_matches_matrix_var<true>();
  normal_lpdf_test::expect_varmat_matches_matrix_var<false>();
}
