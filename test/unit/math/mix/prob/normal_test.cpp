#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/util.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

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

template <typename T, stan::require_std_vector_vt<stan::is_var, T>* = nullptr>
void append_adjoints(std::vector<double>& out, const T& x) {
  for (const auto& e : x) {
    out.push_back(e.adj());
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

// Values and gradients from mpmath (40 digits) at the double inputs of
// check_reference.
constexpr double log_sqrt_two_pi = 0.91893853320467274178;
const std::vector<double> d_mu_scalar_sigma{
    -0.1384083044982699, 0.55363321799307961, 0.17301038062283738,
    -0.24221453287197232};
const std::vector<double> d_mu_vector_sigma{
    -1.6, 0.94674556213017744, 0.10330578512396693, -0.86419753086419743};
const std::vector<double> d_sigma_vector_sigma{
    -0.72000000000000004, 0.3959945380063722, -0.43106686701728021,
    -0.43895747599451314};

// The references are correctly rounded; the rounding error of the function
// is a few ulp, and the tolerance is about 70 ulp of the value.
inline void expect_near_rel(const std::vector<double>& x,
                            const std::vector<double>& ref) {
  ASSERT_EQ(x.size(), ref.size());
  for (size_t i = 0; i < ref.size(); ++i) {
    EXPECT_NEAR(x[i], ref[i], 1e-14 * std::max(1.0, std::fabs(ref[i])))
        << "entry " << i;
  }
}

// Concatenate value, y partials, mu partials and sigma partials.
inline std::vector<double> concat(
    double value, std::initializer_list<std::vector<double>> v) {
  std::vector<double> out{value};
  for (const auto& x : v) {
    out.insert(out.end(), x.begin(), x.end());
  }
  return out;
}

inline std::vector<double> negate(std::vector<double> x) {
  for (auto& e : x) {
    e = -e;
  }
  return x;
}

template <bool propto, typename ToAD>
void check_reference(const ToAD& to_ad) {
  Eigen::VectorXd y(4);
  y << -0.3, 1.2, 2.5, 0.0;
  Eigen::VectorXd mu(4);
  mu << 0.1, -0.4, 2.0, 0.7;
  Eigen::VectorXd sigma(4);
  sigma << 0.5, 1.3, 2.2, 0.9;
  const double add = propto ? 4 * log_sqrt_two_pi : 0.0;
  const double value_scalar_sigma = -6.3968830540223898 + add;
  const double value_vector_sigma = -5.3337600932202835 + add;
  const double d_sigma_scalar_sigma = -1.6486871565235091;
  auto data_y = [&]() { return y; };
  auto ad_sigma = []() { return stan::math::var(1.7); };

  // autodiff mu vector, scalar sigma
  expect_near_rel(
      value_and_adjoints<propto>(
          data_y, [&]() { return to_ad(mu); }, ad_sigma),
      concat(value_scalar_sigma, {d_mu_scalar_sigma, {d_sigma_scalar_sigma}}));
  // autodiff y vector, data mu vector, scalar sigma
  expect_near_rel(value_and_adjoints<propto>([&]() { return to_ad(y); },
                                             [&]() { return mu; }, ad_sigma),
                  concat(value_scalar_sigma,
                         {negate(d_mu_scalar_sigma), {d_sigma_scalar_sigma}}));
  // autodiff mu and sigma vectors
  expect_near_rel(
      value_and_adjoints<propto>(
          data_y, [&]() { return to_ad(mu); }, [&]() { return to_ad(sigma); }),
      concat(value_vector_sigma, {d_mu_vector_sigma, d_sigma_vector_sigma}));
  // all autodiff vectors
  expect_near_rel(
      value_and_adjoints<propto>([&]() { return to_ad(y); },
                                 [&]() { return to_ad(mu); },
                                 [&]() { return to_ad(sigma); }),
      concat(value_vector_sigma, {negate(d_mu_vector_sigma), d_mu_vector_sigma,
                                  d_sigma_vector_sigma}));
}

template <bool propto>
void check_reference_all_types() {
  using stan::math::var;
  check_reference<propto>(
      [](const Eigen::VectorXd& x) { return stan::math::to_var(x); });
  check_reference<propto>([](const Eigen::VectorXd& x) {
    return stan::math::var_value<Eigen::VectorXd>(x);
  });
  check_reference<propto>([](const Eigen::VectorXd& x) {
    return stan::math::var_value<Eigen::RowVectorXd>(x.transpose());
  });
  check_reference<propto>([](const Eigen::VectorXd& x) {
    Eigen::Matrix<var, 1, -1> r = stan::math::to_var(x).transpose();
    return r;
  });
  check_reference<propto>([](const Eigen::VectorXd& x) {
    return std::vector<var>(x.data(), x.data() + x.size());
  });
}
}  // namespace normal_lpdf_test

TEST_F(AgradRev, mathMixScalFun_normal_lpdf_reference) {
  normal_lpdf_test::check_reference_all_types<true>();
  normal_lpdf_test::check_reference_all_types<false>();
}

TEST_F(AgradRev, mathMixScalFun_normal_lpdf_messages) {
  using stan::math::normal_lpdf;
  using stan::math::to_var;
  using stan::math::var;
  using fd = stan::math::fvar<double>;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  Eigen::VectorXd y(3);
  y << 0.5, -0.2, 1.0;
  Eigen::VectorXd mu(3);
  mu << 0.1, 0.0, -0.3;
  Eigen::VectorXd sigma(3);
  sigma << 1.0, 0.7, 1.5;
  Eigen::VectorXd y_nan = y;
  y_nan[1] = nan;
  Eigen::VectorXd mu_inf = mu;
  mu_inf[0] = -inf;
  Eigen::VectorXd sigma_neg = sigma;
  sigma_neg[2] = -1.0;
  const char* y_msg
      = "normal_lpdf: Random variable[2] is nan, but must be not nan!";
  const char* mu_msg
      = "normal_lpdf: Location parameter[1] is -inf, but must be finite!";
  const char* sigma_msg
      = "normal_lpdf: Scale parameter[3] is -1, but must be positive!";
  using vm = stan::math::var_value<Eigen::VectorXd>;
  using vs = std::vector<var>;
  auto to_vs = [](const Eigen::VectorXd& x) {
    return vs(x.data(), x.data() + x.size());
  };

  // the partials of mu hold the scaled residuals
  EXPECT_THROW_MSG(normal_lpdf(y_nan, to_var(mu), var(1.0)), std::domain_error,
                   y_msg);
  EXPECT_THROW_MSG(normal_lpdf(y_nan, vm(mu), var(1.0)), std::domain_error,
                   y_msg);
  EXPECT_THROW_MSG(normal_lpdf(y, vm(mu_inf), var(1.0)), std::domain_error,
                   mu_msg);
  // the partials of y hold the scaled residuals
  EXPECT_THROW_MSG(normal_lpdf(to_var(y), mu_inf, var(1.0)), std::domain_error,
                   mu_msg);
  EXPECT_THROW_MSG(normal_lpdf(to_vs(y_nan), mu, 2.0), std::domain_error,
                   y_msg);
  // autodiff vector sigma (the loop that writes two partials)
  EXPECT_THROW_MSG(normal_lpdf(y_nan, to_vs(mu), to_vs(sigma)),
                   std::domain_error, y_msg);
  EXPECT_THROW_MSG(normal_lpdf(to_var(y), to_var(mu_inf), to_var(sigma)),
                   std::domain_error, mu_msg);
  EXPECT_THROW_MSG(normal_lpdf(to_var(y), to_var(mu), to_var(sigma_neg)),
                   std::domain_error, sigma_msg);
  // no partial of size N other than that of sigma
  EXPECT_THROW_MSG(normal_lpdf(y_nan, mu, to_var(sigma)), std::domain_error,
                   y_msg);
  EXPECT_THROW_MSG(normal_lpdf(y, mu_inf, var(1.0)), std::domain_error, mu_msg);
  // an invalid scalar sigma
  EXPECT_THROW_MSG(normal_lpdf(y, vm(mu), var(-1.0)), std::domain_error,
                   "normal_lpdf: Scale parameter is -1, but must be positive!");
  // forward mode
  Eigen::Matrix<fd, -1, 1> mu_fd = mu.cast<fd>();
  EXPECT_THROW_MSG(normal_lpdf(y_nan, mu_fd, fd(1.0)), std::domain_error,
                   y_msg);
  Eigen::Matrix<fd, -1, 1> sigma_fd = sigma.cast<fd>();
  EXPECT_THROW_MSG(normal_lpdf(y, mu_inf, sigma_fd), std::domain_error, mu_msg);
  Eigen::Matrix<fd, -1, 1> sigma_neg_fd = sigma_neg.cast<fd>();
  EXPECT_THROW_MSG(normal_lpdf(y, mu, sigma_neg_fd), std::domain_error,
                   sigma_msg);
}

TEST_F(AgradRev, mathMixScalFun_normal_lpdf_infinite_y) {
  using stan::math::var;
  const double inf = std::numeric_limits<double>::infinity();
  Eigen::VectorXd y(3);
  y << 0.5, inf, 1.0;
  Eigen::VectorXd mu = Eigen::VectorXd::Zero(3);
  var lp = stan::math::normal_lpdf(y, stan::math::to_var(mu), var(1.0));
  EXPECT_EQ(-inf, lp.val());
  lp = stan::math::normal_lpdf(stan::math::to_var(y), mu, var(1.0));
  EXPECT_EQ(-inf, lp.val());
}
