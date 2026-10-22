#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

// In a vector argument of a discrete cdf, lcdf or lccdf, an element with
// probability 1 contributes nothing, and an element with probability 0 makes
// the result 0 (cdf) or -inf (lcdf, lccdf) with zero derivatives. The tests
// compare f({ordinary, edge}) and f({edge, ordinary}) with f({ordinary}) in
// reverse and forward mode.

namespace {

using stan::math::fvar;
using stan::math::var;

struct value_and_grad {
  double value;
  std::vector<double> grad;
};

template <typename F>
value_and_grad eval_rev(const F& f, const std::vector<int>& n,
                        const std::vector<double>& p) {
  std::vector<var> pv(p.begin(), p.end());
  var y = f(n, pv);
  value_and_grad res{y.val(), {}};
  y.grad();
  for (const auto& x : pv) {
    res.grad.push_back(x.adj());
  }
  stan::math::recover_memory();
  return res;
}

template <typename F>
value_and_grad eval_fwd(const F& f, const std::vector<int>& n,
                        const std::vector<double>& p) {
  value_and_grad res{0, {}};
  for (size_t k = 0; k < p.size(); ++k) {
    std::vector<fvar<double>> pf;
    for (size_t j = 0; j < p.size(); ++j) {
      pf.emplace_back(p[j], j == k ? 1.0 : 0.0);
    }
    fvar<double> y = f(n, pf);
    res.value = y.val();
    res.grad.push_back(y.d_);
  }
  return res;
}

void expect_same(double a, double b, const char* what) {
  if (std::isinf(b)) {
    EXPECT_EQ(a, b) << what;
  } else {
    EXPECT_NEAR(a, b, 1e-12 * std::fmax(1.0, std::fabs(b))) << what;
  }
}

template <typename F>
void expect_probability_one(const F& f, int ordinary, int edge,
                            const std::vector<double>& p) {
  for (int mode = 0; mode < 2; ++mode) {
    auto eval = [&](const std::vector<int>& n) {
      return mode == 0 ? eval_rev(f, n, p) : eval_fwd(f, n, p);
    };
    const value_and_grad alone = eval({ordinary});
    for (const auto& n :
         {std::vector<int>{ordinary, edge}, std::vector<int>{edge, ordinary}}) {
      const value_and_grad both = eval(n);
      expect_same(both.value, alone.value,
                  mode == 0 ? "rev value" : "fwd value");
      for (size_t k = 0; k < p.size(); ++k) {
        expect_same(both.grad[k], alone.grad[k],
                    mode == 0 ? "rev gradient" : "fwd gradient");
      }
    }
  }
}

template <typename F>
void expect_probability_zero(const F& f, int ordinary, int edge,
                             const std::vector<double>& p, double result) {
  for (int mode = 0; mode < 2; ++mode) {
    for (const auto& n :
         {std::vector<int>{ordinary, edge}, std::vector<int>{edge, ordinary}}) {
      const value_and_grad both
          = mode == 0 ? eval_rev(f, n, p) : eval_fwd(f, n, p);
      EXPECT_EQ(both.value, result) << (mode == 0 ? "rev" : "fwd");
      for (size_t k = 0; k < p.size(); ++k) {
        EXPECT_EQ(both.grad[k], 0.0) << (mode == 0 ? "rev" : "fwd");
      }
    }
  }
}

constexpr double NEG_INF = -std::numeric_limits<double>::infinity();
constexpr int INT_MAX_N = std::numeric_limits<int>::max();

}  // namespace

TEST(ProbDiscreteCdfsEdgeElements, bernoulli) {
  auto cdf = [](const auto& n, const auto& p) {
    return stan::math::bernoulli_cdf(n, p[0]);
  };
  auto lcdf = [](const auto& n, const auto& p) {
    return stan::math::bernoulli_lcdf(n, p[0]);
  };
  auto lccdf = [](const auto& n, const auto& p) {
    return stan::math::bernoulli_lccdf(n, p[0]);
  };
  const std::vector<double> p{0.3};
  expect_probability_one(cdf, 0, 1, p);
  expect_probability_one(lcdf, 0, 1, p);
  expect_probability_one(lccdf, 0, -1, p);
  expect_probability_zero(cdf, 0, -1, p, 0.0);
  expect_probability_zero(lcdf, 0, -1, p, NEG_INF);
  expect_probability_zero(lccdf, 0, 1, p, NEG_INF);
}

TEST(ProbDiscreteCdfsEdgeElements, binomial) {
  auto cdf = [](const auto& n, const auto& p) {
    return stan::math::binomial_cdf(n, 10, p[0]);
  };
  auto lcdf = [](const auto& n, const auto& p) {
    return stan::math::binomial_lcdf(n, 10, p[0]);
  };
  auto lccdf = [](const auto& n, const auto& p) {
    return stan::math::binomial_lccdf(n, 10, p[0]);
  };
  const std::vector<double> p{0.3};
  expect_probability_one(cdf, 3, 10, p);
  expect_probability_one(lcdf, 3, 10, p);
  expect_probability_one(lccdf, 3, -1, p);
  expect_probability_zero(cdf, 3, -1, p, 0.0);
  expect_probability_zero(lcdf, 3, -1, p, NEG_INF);
  expect_probability_zero(lccdf, 3, 10, p, NEG_INF);
}

TEST(ProbDiscreteCdfsEdgeElements, beta_binomial) {
  auto cdf = [](const auto& n, const auto& p) {
    return stan::math::beta_binomial_cdf(n, 10, p[0], p[1]);
  };
  auto lcdf = [](const auto& n, const auto& p) {
    return stan::math::beta_binomial_lcdf(n, 10, p[0], p[1]);
  };
  auto lccdf = [](const auto& n, const auto& p) {
    return stan::math::beta_binomial_lccdf(n, 10, p[0], p[1]);
  };
  const std::vector<double> p{2.5, 1.5};
  expect_probability_one(cdf, 3, 10, p);
  expect_probability_one(lcdf, 3, 10, p);
  expect_probability_one(lccdf, 3, -1, p);
  expect_probability_zero(cdf, 3, -1, p, 0.0);
  expect_probability_zero(lcdf, 3, -1, p, NEG_INF);
  expect_probability_zero(lccdf, 3, 10, p, NEG_INF);
}

TEST(ProbDiscreteCdfsEdgeElements, neg_binomial) {
  auto cdf = [](const auto& n, const auto& p) {
    return stan::math::neg_binomial_cdf(n, p[0], p[1]);
  };
  auto lcdf = [](const auto& n, const auto& p) {
    return stan::math::neg_binomial_lcdf(n, p[0], p[1]);
  };
  auto lccdf = [](const auto& n, const auto& p) {
    return stan::math::neg_binomial_lccdf(n, p[0], p[1]);
  };
  const std::vector<double> p{2.5, 1.5};
  expect_probability_one(cdf, 3, INT_MAX_N, p);
  expect_probability_one(lcdf, 3, INT_MAX_N, p);
  expect_probability_one(lccdf, 3, -1, p);
  expect_probability_zero(cdf, 3, -1, p, 0.0);
  expect_probability_zero(lcdf, 3, -1, p, NEG_INF);
  expect_probability_zero(lccdf, 3, INT_MAX_N, p, NEG_INF);
}

TEST(ProbDiscreteCdfsEdgeElements, neg_binomial_2) {
  auto cdf = [](const auto& n, const auto& p) {
    return stan::math::neg_binomial_2_cdf(n, p[0], p[1]);
  };
  auto lcdf = [](const auto& n, const auto& p) {
    return stan::math::neg_binomial_2_lcdf(n, p[0], p[1]);
  };
  auto lccdf = [](const auto& n, const auto& p) {
    return stan::math::neg_binomial_2_lccdf(n, p[0], p[1]);
  };
  const std::vector<double> p{4.0, 2.0};
  expect_probability_one(cdf, 3, INT_MAX_N, p);
  expect_probability_one(lcdf, 3, INT_MAX_N, p);
  expect_probability_one(lccdf, 3, -1, p);
  expect_probability_zero(cdf, 3, -1, p, 0.0);
  expect_probability_zero(lcdf, 3, -1, p, NEG_INF);
  expect_probability_zero(lccdf, 3, INT_MAX_N, p, NEG_INF);
}

TEST(ProbDiscreteCdfsEdgeElements, beta_neg_binomial) {
  auto cdf = [](const auto& n, const auto& p) {
    return stan::math::beta_neg_binomial_cdf(n, p[0], p[1], p[2]);
  };
  auto lcdf = [](const auto& n, const auto& p) {
    return stan::math::beta_neg_binomial_lcdf(n, p[0], p[1], p[2]);
  };
  auto lccdf = [](const auto& n, const auto& p) {
    return stan::math::beta_neg_binomial_lccdf(n, p[0], p[1], p[2]);
  };
  const std::vector<double> p{3.0, 2.5, 1.5};
  expect_probability_one(cdf, 3, INT_MAX_N, p);
  expect_probability_one(lcdf, 3, INT_MAX_N, p);
  expect_probability_one(lccdf, 3, -1, p);
  expect_probability_zero(cdf, 3, -1, p, 0.0);
  expect_probability_zero(lcdf, 3, -1, p, NEG_INF);
  expect_probability_zero(lccdf, 3, INT_MAX_N, p, NEG_INF);
}

TEST(ProbDiscreteCdfsEdgeElements, yule_simon) {
  auto cdf = [](const auto& n, const auto& p) {
    return stan::math::yule_simon_cdf(n, p[0]);
  };
  auto lcdf = [](const auto& n, const auto& p) {
    return stan::math::yule_simon_lcdf(n, p[0]);
  };
  auto lccdf = [](const auto& n, const auto& p) {
    return stan::math::yule_simon_lccdf(n, p[0]);
  };
  const std::vector<double> p{1.5};
  expect_probability_one(cdf, 3, INT_MAX_N, p);
  expect_probability_one(lcdf, 3, INT_MAX_N, p);
  expect_probability_one(lccdf, 3, 0, p);
  expect_probability_zero(cdf, 3, 0, p, 0.0);
  expect_probability_zero(lcdf, 3, 0, p, NEG_INF);
  expect_probability_zero(lccdf, 3, INT_MAX_N, p, NEG_INF);
}
