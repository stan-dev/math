#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <functional>
#include <limits>

namespace {

struct grad_case {
  double y;
  double mu;
  double lambda;
  double d_y;
  double d_mu;
  double d_lambda;
};

template <typename F>
void check_grads(const F& f, const grad_case& c, double tol) {
  stan::math::var y = c.y;
  stan::math::var mu = c.mu;
  stan::math::var lambda = c.lambda;
  stan::math::var out = f(y, mu, lambda);
  out.grad();
  EXPECT_NEAR(c.d_y, y.adj(), tol);
  EXPECT_NEAR(c.d_mu, mu.adj(), tol);
  EXPECT_NEAR(c.d_lambda, lambda.adj(), tol);
  stan::math::recover_memory();
}

}  // namespace

// reference values from mpmath at 60 digits
TEST_F(AgradRev, inv_gaussian_lpdf_gradients) {
  auto f = [](const auto& y, const auto& mu, const auto& lambda) {
    return stan::math::inv_gaussian_lpdf(y, mu, lambda);
  };
  check_grads(f,
              {1.2, 0.5, 2.0, -4.5555555555555556, 11.199999999999999,
               -0.56666666666666659},
              1e-12);
  check_grads(f,
              {0.3, 1.0, 5.0, 20.27777777777778, -3.5000000000000001,
               -0.71666666666666672},
              1e-12);
}

TEST_F(AgradRev, inv_gaussian_lcdf_gradients) {
  auto f = [](const auto& y, const auto& mu, const auto& lambda) {
    return stan::math::inv_gaussian_lcdf(y, mu, lambda);
  };
  check_grads(f,
              {1.2, 0.5, 2.0, 0.085383591493919688, -0.27616885497947743,
               0.017812058848517546},
              1e-12);
  check_grads(f,
              {0.3, 1.0, 5.0, 27.233870297028338, -3.6491778320016701,
               -0.90419665142136623},
              1e-11);
}

TEST_F(AgradRev, inv_gaussian_lccdf_gradients) {
  auto f = [](const auto& y, const auto& mu, const auto& lambda) {
    return stan::math::inv_gaussian_lccdf(y, mu, lambda);
  };
  check_grads(f,
              {1.2, 0.5, 2.0, -4.5530760732154196, 14.72669143770909,
               -0.94982721549802085},
              1e-11);
  check_grads(f,
              {0.3, 1.0, 5.0, -0.09179211684684445, 0.012299601720089049,
               0.003047606666792857},
              1e-12);
}

// A batch of N observations produces one node on the autodiff tape, so tape
// growth beyond the input vars must not depend on N.
TEST_F(AgradRev, inv_gaussian_one_node_per_call) {
  using stan::math::var;
  using stan::math::vector_v;

  auto tape_growth = [](int n, auto&& f) {
    vector_v y = Eigen::VectorXd::LinSpaced(n, 0.4, 2.0);
    vector_v mu = Eigen::VectorXd::Constant(n, 1.0);
    vector_v lambda = Eigen::VectorXd::Constant(n, 2.0);
    std::size_t before
        = stan::math::ChainableStack::instance_->var_stack_.size();
    var out = f(y, mu, lambda);
    std::size_t after
        = stan::math::ChainableStack::instance_->var_stack_.size();
    return after - before;
  };

  auto lpdf = [](const auto& y, const auto& mu, const auto& lambda) {
    return stan::math::inv_gaussian_lpdf(y, mu, lambda);
  };
  auto lcdf = [](const auto& y, const auto& mu, const auto& lambda) {
    return stan::math::inv_gaussian_lcdf(y, mu, lambda);
  };
  auto lccdf = [](const auto& y, const auto& mu, const auto& lambda) {
    return stan::math::inv_gaussian_lccdf(y, mu, lambda);
  };
  auto cdf = [](const auto& y, const auto& mu, const auto& lambda) {
    return stan::math::inv_gaussian_cdf(y, mu, lambda);
  };

  for (auto&& f :
       {std::function<var(const vector_v&, const vector_v&, const vector_v&)>(
            lpdf),
        std::function<var(const vector_v&, const vector_v&, const vector_v&)>(
            lcdf),
        std::function<var(const vector_v&, const vector_v&, const vector_v&)>(
            lccdf),
        std::function<var(const vector_v&, const vector_v&, const vector_v&)>(
            cdf)}) {
    std::size_t g10 = tape_growth(10, f);
    std::size_t g1000 = tape_growth(1000, f);
    EXPECT_EQ(g10, g1000);
    stan::math::recover_memory();
  }
}

TEST_F(AgradRev, inv_gaussian_vectorized_matches_scalar) {
  using stan::math::var;
  using stan::math::vector_v;
  Eigen::VectorXd y_d(4);
  y_d << 0.3, 0.8, 1.5, 3.0;
  Eigen::VectorXd mu_d(4);
  mu_d << 0.5, 1.0, 1.0, 2.0;
  Eigen::VectorXd lambda_d(4);
  lambda_d << 2.0, 5.0, 1.0, 0.7;

  vector_v y = y_d;
  vector_v mu = mu_d;
  vector_v lambda = lambda_d;
  var out = stan::math::inv_gaussian_lccdf(y, mu, lambda);
  out.grad();
  Eigen::VectorXd vec_adj_y(4);
  for (int i = 0; i < 4; ++i) {
    vec_adj_y(i) = y(i).adj();
  }
  double vec_val = out.val();
  stan::math::recover_memory();

  double scalar_sum = 0;
  Eigen::VectorXd scalar_adj_y(4);
  for (int i = 0; i < 4; ++i) {
    var yi = y_d(i);
    var out_i = stan::math::inv_gaussian_lccdf(yi, mu_d(i), lambda_d(i));
    out_i.grad();
    scalar_sum += out_i.val();
    scalar_adj_y(i) = yi.adj();
    stan::math::recover_memory();
  }

  EXPECT_FLOAT_EQ(scalar_sum, vec_val);
  for (int i = 0; i < 4; ++i) {
    EXPECT_FLOAT_EQ(scalar_adj_y(i), vec_adj_y(i));
  }
}

// Boundary adjoints are asserted directly here; expect_ad's finite-difference
// stencil steps off the support at y == 0 and is NaN at y == inf.
TEST_F(AgradRev, inv_gaussian_boundary_partials_are_zero) {
  using stan::math::var;
  double inf = std::numeric_limits<double>::infinity();

  auto check = [](auto f, double yv, const char* what) {
    var y = yv;
    var mu = 1.0;
    var lambda = 2.0;
    var out = f(y, mu, lambda);
    out.grad();
    EXPECT_FLOAT_EQ(0.0, y.adj()) << what << " d/dy at y=" << yv;
    EXPECT_FLOAT_EQ(0.0, mu.adj()) << what << " d/dmu at y=" << yv;
    EXPECT_FLOAT_EQ(0.0, lambda.adj()) << what << " d/dlambda at y=" << yv;
    stan::math::recover_memory();
  };

  auto lpdf = [](const auto& y, const auto& mu, const auto& lambda) {
    return stan::math::inv_gaussian_lpdf(y, mu, lambda);
  };
  auto cdf = [](const auto& y, const auto& mu, const auto& lambda) {
    return stan::math::inv_gaussian_cdf(y, mu, lambda);
  };
  auto lcdf = [](const auto& y, const auto& mu, const auto& lambda) {
    return stan::math::inv_gaussian_lcdf(y, mu, lambda);
  };
  auto lccdf = [](const auto& y, const auto& mu, const auto& lambda) {
    return stan::math::inv_gaussian_lccdf(y, mu, lambda);
  };

  for (double b : {0.0, inf}) {
    check(lpdf, b, "lpdf");
    check(cdf, b, "cdf");
    check(lcdf, b, "lcdf");
    check(lccdf, b, "lccdf");
  }
}

// The log probability saturates to -inf at an interior y in two ways: the
// guarded survivor difference underflows deep in the upper tail, and
// lambda / y overflows for representable arguments. Both trigger the same
// mask, which zeroes every partial.
TEST_F(AgradRev, inv_gaussian_saturated_partials_are_zero) {
  using stan::math::var;
  double inf = std::numeric_limits<double>::infinity();

  auto check = [](auto f, double y_dbl, double mu_dbl, double lambda_dbl,
                  double val, const char* what) {
    var y = y_dbl;
    var mu = mu_dbl;
    var lambda = lambda_dbl;
    var out = f(y, mu, lambda);
    EXPECT_FLOAT_EQ(val, out.val()) << what;
    out.grad();
    EXPECT_FLOAT_EQ(0.0, y.adj()) << what;
    EXPECT_FLOAT_EQ(0.0, mu.adj()) << what;
    EXPECT_FLOAT_EQ(0.0, lambda.adj()) << what;
    stan::math::recover_memory();
  };

  auto lcdf = [](const auto& y, const auto& mu, const auto& lambda) {
    return stan::math::inv_gaussian_lcdf(y, mu, lambda);
  };
  auto cdf = [](const auto& y, const auto& mu, const auto& lambda) {
    return stan::math::inv_gaussian_cdf(y, mu, lambda);
  };
  auto lccdf = [](const auto& y, const auto& mu, const auto& lambda) {
    return stan::math::inv_gaussian_lccdf(y, mu, lambda);
  };

  // lambda / y = 1e310 overflows; y is far below mu, so the lower tail
  // saturates
  check(lcdf, 1e-10, 1.0, 1e300, -inf, "lcdf");
  check(cdf, 1e-10, 1.0, 1e300, 0.0, "cdf");
  // y is above mu, so the survivor saturates
  check(lccdf, 2e-3, 1e-3, 1e308, -inf, "lccdf");
  // survivor underflow: the two lccdf terms agree to the last bit
  check(lccdf, 1e-2, 1e-3, 1e14, -inf, "lccdf underflow");
  check(lccdf, 1.0, 1e-3, 1e14, -inf, "lccdf underflow");
  check(lccdf, 10.0, 1e-3, 1e14, -inf, "lccdf underflow");
}

// A boundary element leaves the other elements' adjoints alone; the
// container values themselves are pinned in the prim suite.
TEST_F(AgradRev, inv_gaussian_boundary_adjoints_are_elementwise) {
  using stan::math::var;
  using stan::math::vector_v;
  double inf = std::numeric_limits<double>::infinity();

  auto check = [](auto f, double boundary, double finite) {
    vector_v y_v(2);
    y_v << boundary, finite;
    var mu = 1.0;
    var lambda = 2.0;
    var out = f(y_v, mu, lambda);
    out.grad();
    double adj_boundary = y_v(0).adj();
    double adj_finite = y_v(1).adj();
    double adj_mu = mu.adj();
    double adj_lambda = lambda.adj();
    stan::math::recover_memory();

    var y_s = finite;
    var mu_s = 1.0;
    var lambda_s = 2.0;
    var out_s = f(y_s, mu_s, lambda_s);
    out_s.grad();
    EXPECT_FLOAT_EQ(0.0, adj_boundary);
    EXPECT_FLOAT_EQ(y_s.adj(), adj_finite);
    EXPECT_FLOAT_EQ(mu_s.adj(), adj_mu);
    EXPECT_FLOAT_EQ(lambda_s.adj(), adj_lambda);
    stan::math::recover_memory();
  };

  check(
      [](const auto& y, const auto& mu, const auto& lambda) {
        return stan::math::inv_gaussian_lccdf(y, mu, lambda);
      },
      0.0, 5.0);
  check(
      [](const auto& y, const auto& mu, const auto& lambda) {
        return stan::math::inv_gaussian_lcdf(y, mu, lambda);
      },
      inf, 0.5);
}
