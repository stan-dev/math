#include <test/unit/math/test_ad.hpp>
#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <type_traits>
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

// Parameter j of k, with p holding the k parameters one after the other, each
// as a vector of p.size() / k elements
template <typename T>
std::vector<T> part(const std::vector<T>& p, size_t k, size_t j) {
  const size_t m = p.size() / k;
  return std::vector<T>(p.begin() + j * m, p.begin() + (j + 1) * m);
}

// Vector parameters: q_ord and q_edge are the parameters of the ordinary and
// of the edge element. The derivatives with respect to the parameters of the
// edge element must be exactly 0.
template <typename F>
void expect_probability_one_vec(const F& f, int ordinary, int edge,
                                const std::vector<double>& q_ord,
                                const std::vector<double>& q_edge) {
  const size_t k = q_ord.size();
  for (int mode = 0; mode < 2; ++mode) {
    auto eval = [&](const std::vector<int>& n, const std::vector<double>& p) {
      return mode == 0 ? eval_rev(f, n, p) : eval_fwd(f, n, p);
    };
    const value_and_grad alone = eval({ordinary}, q_ord);
    for (size_t pos_ord = 0; pos_ord < 2; ++pos_ord) {
      const size_t pos_edge = 1 - pos_ord;
      std::vector<int> n(2);
      n[pos_ord] = ordinary;
      n[pos_edge] = edge;
      std::vector<double> p(2 * k);
      for (size_t j = 0; j < k; ++j) {
        p[2 * j + pos_ord] = q_ord[j];
        p[2 * j + pos_edge] = q_edge[j];
      }
      const value_and_grad both = eval(n, p);
      const char* what = mode == 0 ? "rev" : "fwd";
      expect_same(both.value, alone.value, what);
      for (size_t j = 0; j < k; ++j) {
        expect_same(both.grad[2 * j + pos_ord], alone.grad[j], what);
        EXPECT_EQ(both.grad[2 * j + pos_edge], 0.0) << what;
      }
    }
  }
}

// A scalar edge count with vector parameters: the result is the given value
// with zero derivatives
template <typename F>
void expect_scalar_edge(const F& f, int edge, const std::vector<double>& p,
                        double result) {
  auto g = [&](const std::vector<int>& n, const auto& q) { return f(n[0], q); };
  for (int mode = 0; mode < 2; ++mode) {
    const value_and_grad res
        = mode == 0 ? eval_rev(g, {edge}, p) : eval_fwd(g, {edge}, p);
    EXPECT_EQ(res.value, result) << (mode == 0 ? "rev" : "fwd");
    for (size_t k = 0; k < p.size(); ++k) {
      EXPECT_EQ(res.grad[k], 0.0) << (mode == 0 ? "rev" : "fwd");
    }
  }
}

// poisson_binomial: element e has the k probabilities p[j * m + e],
// j = 0, ..., k - 1, with m = p.size() / k elements
template <typename T>
std::vector<Eigen::Matrix<T, -1, 1>> make_theta(const std::vector<T>& p,
                                                size_t k) {
  const size_t m = p.size() / k;
  std::vector<Eigen::Matrix<T, -1, 1>> theta;
  for (size_t e = 0; e < m; ++e) {
    Eigen::Matrix<T, -1, 1> th(k);
    for (size_t j = 0; j < k; ++j) {
      th(j) = p[j * m + e];
    }
    theta.push_back(th);
  }
  return theta;
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

TEST(ProbDiscreteCdfsEdgeElements, poisson_binomial) {
  auto theta = [](const auto& y, const auto& p) {
    using T = std::decay_t<decltype(p[0])>;
    Eigen::Matrix<T, -1, 1> th(3);
    th << p[0], p[1], p[2];
    return std::vector<Eigen::Matrix<T, -1, 1>>(y.size(), th);
  };
  auto cdf = [&](const auto& y, const auto& p) {
    return stan::math::poisson_binomial_cdf(y, theta(y, p));
  };
  auto lcdf = [&](const auto& y, const auto& p) {
    return stan::math::poisson_binomial_lcdf(y, theta(y, p));
  };
  auto lccdf = [&](const auto& y, const auto& p) {
    return stan::math::poisson_binomial_lccdf(y, theta(y, p));
  };
  const std::vector<double> p{0.3, 0.6, 0.8};
  expect_probability_one(cdf, 1, 3, p);
  expect_probability_one(lcdf, 1, 3, p);
  expect_probability_zero(lccdf, 1, 3, p, NEG_INF);
}

TEST(ProbDiscreteCdfsEdgeElements, vector_parameters) {
  using stan::math::bernoulli_cdf;
  using stan::math::bernoulli_lccdf;
  using stan::math::bernoulli_lcdf;
  const std::vector<double> t_ord{0.3}, t_edge{0.6};
  auto b_cdf = [](const auto& n, const auto& p) { return bernoulli_cdf(n, p); };
  auto b_lcdf
      = [](const auto& n, const auto& p) { return bernoulli_lcdf(n, p); };
  auto b_lccdf
      = [](const auto& n, const auto& p) { return bernoulli_lccdf(n, p); };
  expect_probability_one_vec(b_cdf, 0, 1, t_ord, t_edge);
  expect_probability_one_vec(b_lcdf, 0, 1, t_ord, t_edge);
  expect_probability_one_vec(b_lccdf, 0, -1, t_ord, t_edge);
  expect_probability_zero(b_lccdf, 0, 1, {0.3, 0.6}, NEG_INF);

  auto bin_lccdf = [](const auto& n, const auto& p) {
    return stan::math::binomial_lccdf(n, 10, p);
  };
  expect_probability_one_vec(bin_lccdf, 3, -1, t_ord, t_edge);
  expect_probability_zero(bin_lccdf, 3, 10, {0.3, 0.6}, NEG_INF);

  const std::vector<double> ab_ord{2.5, 1.5}, ab_edge{3.5, 0.8};
  const std::vector<double> ab_both{2.5, 3.5, 1.5, 0.8};
  auto bb_lccdf = [](const auto& n, const auto& p) {
    return stan::math::beta_binomial_lccdf(n, 10, part(p, 2, 0), part(p, 2, 1));
  };
  expect_probability_one_vec(bb_lccdf, 3, -1, ab_ord, ab_edge);
  expect_probability_zero(bb_lccdf, 3, 10, ab_both, NEG_INF);

  auto nb_cdf = [](const auto& n, const auto& p) {
    return stan::math::neg_binomial_cdf(n, part(p, 2, 0), part(p, 2, 1));
  };
  auto nb_lcdf = [](const auto& n, const auto& p) {
    return stan::math::neg_binomial_lcdf(n, part(p, 2, 0), part(p, 2, 1));
  };
  auto nb_lccdf = [](const auto& n, const auto& p) {
    return stan::math::neg_binomial_lccdf(n, part(p, 2, 0), part(p, 2, 1));
  };
  expect_probability_one_vec(nb_cdf, 3, INT_MAX_N, ab_ord, ab_edge);
  expect_probability_one_vec(nb_lcdf, 3, INT_MAX_N, ab_ord, ab_edge);
  expect_probability_one_vec(nb_lccdf, 3, -1, ab_ord, ab_edge);
  expect_probability_zero(nb_lccdf, 3, INT_MAX_N, ab_both, NEG_INF);

  const std::vector<double> mp_ord{4.0, 2.0}, mp_edge{1.5, 3.0};
  auto nb2_cdf = [](const auto& n, const auto& p) {
    return stan::math::neg_binomial_2_cdf(n, part(p, 2, 0), part(p, 2, 1));
  };
  auto nb2_lcdf = [](const auto& n, const auto& p) {
    return stan::math::neg_binomial_2_lcdf(n, part(p, 2, 0), part(p, 2, 1));
  };
  auto nb2_lccdf = [](const auto& n, const auto& p) {
    return stan::math::neg_binomial_2_lccdf(n, part(p, 2, 0), part(p, 2, 1));
  };
  expect_probability_one_vec(nb2_cdf, 3, INT_MAX_N, mp_ord, mp_edge);
  expect_probability_one_vec(nb2_lcdf, 3, INT_MAX_N, mp_ord, mp_edge);
  expect_probability_one_vec(nb2_lccdf, 3, -1, mp_ord, mp_edge);

  const std::vector<double> rab_ord{3.0, 2.5, 1.5}, rab_edge{2.0, 3.5, 1.2};
  const std::vector<double> rab_both{3.0, 2.0, 2.5, 3.5, 1.5, 1.2};
  auto bnb_cdf = [](const auto& n, const auto& p) {
    return stan::math::beta_neg_binomial_cdf(n, part(p, 3, 0), part(p, 3, 1),
                                             part(p, 3, 2));
  };
  auto bnb_lcdf = [](const auto& n, const auto& p) {
    return stan::math::beta_neg_binomial_lcdf(n, part(p, 3, 0), part(p, 3, 1),
                                              part(p, 3, 2));
  };
  auto bnb_lccdf = [](const auto& n, const auto& p) {
    return stan::math::beta_neg_binomial_lccdf(n, part(p, 3, 0), part(p, 3, 1),
                                               part(p, 3, 2));
  };
  expect_probability_one_vec(bnb_cdf, 3, INT_MAX_N, rab_ord, rab_edge);
  expect_probability_one_vec(bnb_lcdf, 3, INT_MAX_N, rab_ord, rab_edge);
  expect_probability_one_vec(bnb_lccdf, 3, -1, rab_ord, rab_edge);
  expect_probability_zero(bnb_lccdf, 3, INT_MAX_N, rab_both, NEG_INF);

  const std::vector<double> a_ord{1.5}, a_edge{0.7};
  auto ys_cdf = [](const auto& n, const auto& p) {
    return stan::math::yule_simon_cdf(n, p);
  };
  auto ys_lcdf = [](const auto& n, const auto& p) {
    return stan::math::yule_simon_lcdf(n, p);
  };
  auto ys_lccdf = [](const auto& n, const auto& p) {
    return stan::math::yule_simon_lccdf(n, p);
  };
  expect_probability_one_vec(ys_cdf, 3, INT_MAX_N, a_ord, a_edge);
  expect_probability_one_vec(ys_lcdf, 3, INT_MAX_N, a_ord, a_edge);
  expect_probability_one_vec(ys_lccdf, 3, 0, a_ord, a_edge);

  // A scalar edge count with vector parameters
  expect_scalar_edge(b_lccdf, -1, {0.3, 0.6}, 0.0);
  expect_scalar_edge(bin_lccdf, -1, {0.3, 0.6}, 0.0);
  expect_scalar_edge(bb_lccdf, -1, ab_both, 0.0);
  expect_scalar_edge(nb_cdf, INT_MAX_N, ab_both, 1.0);
  expect_scalar_edge(nb_lcdf, INT_MAX_N, ab_both, 0.0);
  expect_scalar_edge(nb_lccdf, -1, ab_both, 0.0);
  expect_scalar_edge(nb2_cdf, INT_MAX_N, {4.0, 1.5, 2.0, 3.0}, 1.0);
  expect_scalar_edge(nb2_lcdf, INT_MAX_N, {4.0, 1.5, 2.0, 3.0}, 0.0);
  expect_scalar_edge(nb2_lccdf, -1, {4.0, 1.5, 2.0, 3.0}, 0.0);
  expect_scalar_edge(bnb_cdf, INT_MAX_N, rab_both, 1.0);
  expect_scalar_edge(bnb_lcdf, INT_MAX_N, rab_both, 0.0);
  expect_scalar_edge(bnb_lccdf, -1, rab_both, 0.0);
  expect_scalar_edge(ys_cdf, INT_MAX_N, {1.5, 0.7}, 1.0);
  expect_scalar_edge(ys_lcdf, INT_MAX_N, {1.5, 0.7}, 0.0);
  expect_scalar_edge(ys_lccdf, 0, {1.5, 0.7}, 0.0);

  // The binomial and beta_binomial cdf and lcdf with vector parameters
  auto bin_cdf = [](const auto& n, const auto& p) {
    return stan::math::binomial_cdf(n, 10, p);
  };
  auto bin_lcdf = [](const auto& n, const auto& p) {
    return stan::math::binomial_lcdf(n, 10, p);
  };
  auto bb_cdf = [](const auto& n, const auto& p) {
    return stan::math::beta_binomial_cdf(n, 10, part(p, 2, 0), part(p, 2, 1));
  };
  auto bb_lcdf = [](const auto& n, const auto& p) {
    return stan::math::beta_binomial_lcdf(n, 10, part(p, 2, 0), part(p, 2, 1));
  };
  expect_probability_one_vec(bin_cdf, 3, 10, t_ord, t_edge);
  expect_probability_one_vec(bin_lcdf, 3, 10, t_ord, t_edge);
  expect_probability_one_vec(bb_cdf, 3, 10, ab_ord, ab_edge);
  expect_probability_one_vec(bb_lcdf, 3, 10, ab_ord, ab_edge);
  expect_scalar_edge(b_cdf, 1, {0.3, 0.6}, 1.0);
  expect_scalar_edge(b_lcdf, 1, {0.3, 0.6}, 0.0);
  expect_scalar_edge(bin_cdf, 10, {0.3, 0.6}, 1.0);
  expect_scalar_edge(bin_lcdf, 10, {0.3, 0.6}, 0.0);
  expect_scalar_edge(bb_cdf, 10, ab_both, 1.0);
  expect_scalar_edge(bb_lcdf, 10, ab_both, 0.0);

  // Probability 0 with vector parameters
  const std::vector<double> t_both{0.3, 0.6};
  const std::vector<double> mp_both{4.0, 1.5, 2.0, 3.0};
  const std::vector<double> a_both{1.5, 0.7};
  expect_probability_zero(b_cdf, 0, -1, t_both, 0.0);
  expect_probability_zero(b_lcdf, 0, -1, t_both, NEG_INF);
  expect_probability_zero(bin_cdf, 3, -1, t_both, 0.0);
  expect_probability_zero(bin_lcdf, 3, -1, t_both, NEG_INF);
  expect_probability_zero(bb_cdf, 3, -1, ab_both, 0.0);
  expect_probability_zero(bb_lcdf, 3, -1, ab_both, NEG_INF);
  expect_probability_zero(nb_cdf, 3, -1, ab_both, 0.0);
  expect_probability_zero(nb_lcdf, 3, -1, ab_both, NEG_INF);
  expect_probability_zero(nb2_cdf, 3, -1, mp_both, 0.0);
  expect_probability_zero(nb2_lcdf, 3, -1, mp_both, NEG_INF);
  expect_probability_zero(nb2_lccdf, 3, INT_MAX_N, mp_both, NEG_INF);
  expect_probability_zero(bnb_cdf, 3, -1, rab_both, 0.0);
  expect_probability_zero(bnb_lcdf, 3, -1, rab_both, NEG_INF);
  expect_probability_zero(ys_cdf, 3, 0, a_both, 0.0);
  expect_probability_zero(ys_lcdf, 3, 0, a_both, NEG_INF);
  expect_probability_zero(ys_lccdf, 3, INT_MAX_N, a_both, NEG_INF);

  // A scalar edge count with probability 0 and vector parameters
  expect_scalar_edge(b_cdf, -1, t_both, 0.0);
  expect_scalar_edge(b_lcdf, -1, t_both, NEG_INF);
  expect_scalar_edge(b_lccdf, 1, t_both, NEG_INF);
  expect_scalar_edge(bin_cdf, -1, t_both, 0.0);
  expect_scalar_edge(bin_lcdf, -1, t_both, NEG_INF);
  expect_scalar_edge(bin_lccdf, 10, t_both, NEG_INF);
  expect_scalar_edge(bb_cdf, -1, ab_both, 0.0);
  expect_scalar_edge(bb_lcdf, -1, ab_both, NEG_INF);
  expect_scalar_edge(bb_lccdf, 10, ab_both, NEG_INF);
  expect_scalar_edge(nb_cdf, -1, ab_both, 0.0);
  expect_scalar_edge(nb_lcdf, -1, ab_both, NEG_INF);
  expect_scalar_edge(nb_lccdf, INT_MAX_N, ab_both, NEG_INF);
  expect_scalar_edge(nb2_cdf, -1, mp_both, 0.0);
  expect_scalar_edge(nb2_lcdf, -1, mp_both, NEG_INF);
  expect_scalar_edge(nb2_lccdf, INT_MAX_N, mp_both, NEG_INF);
  expect_scalar_edge(bnb_cdf, -1, rab_both, 0.0);
  expect_scalar_edge(bnb_lcdf, -1, rab_both, NEG_INF);
  expect_scalar_edge(bnb_lccdf, INT_MAX_N, rab_both, NEG_INF);
  expect_scalar_edge(ys_cdf, 0, a_both, 0.0);
  expect_scalar_edge(ys_lcdf, 0, a_both, NEG_INF);
  expect_scalar_edge(ys_lccdf, INT_MAX_N, a_both, NEG_INF);
}

TEST(ProbDiscreteCdfsEdgeElements, poisson_binomial_theta_per_element) {
  // Three trials per element
  auto cdf = [](const auto& y, const auto& p) {
    return stan::math::poisson_binomial_cdf(y, make_theta(p, 3));
  };
  auto lcdf = [](const auto& y, const auto& p) {
    return stan::math::poisson_binomial_lcdf(y, make_theta(p, 3));
  };
  auto lccdf = [](const auto& y, const auto& p) {
    return stan::math::poisson_binomial_lccdf(y, make_theta(p, 3));
  };
  const std::vector<double> q_ord{0.3, 0.6, 0.8}, q_edge{0.2, 0.5, 0.9};
  const std::vector<double> q_both{0.3, 0.2, 0.6, 0.5, 0.8, 0.9};
  expect_probability_one_vec(cdf, 1, 3, q_ord, q_edge);
  expect_probability_one_vec(lcdf, 1, 3, q_ord, q_edge);
  expect_probability_zero(lccdf, 1, 3, q_both, NEG_INF);
  expect_scalar_edge(cdf, 3, q_both, 1.0);
  expect_scalar_edge(lcdf, 3, q_both, 0.0);
  expect_scalar_edge(lccdf, 3, q_both, NEG_INF);

  // One trial per element
  auto cdf1 = [](const auto& y, const auto& p) {
    return stan::math::poisson_binomial_cdf(y, make_theta(p, 1));
  };
  auto lcdf1 = [](const auto& y, const auto& p) {
    return stan::math::poisson_binomial_lcdf(y, make_theta(p, 1));
  };
  auto lccdf1 = [](const auto& y, const auto& p) {
    return stan::math::poisson_binomial_lccdf(y, make_theta(p, 1));
  };
  expect_probability_one_vec(cdf1, 0, 1, {0.3}, {0.6});
  expect_probability_one_vec(lcdf1, 0, 1, {0.3}, {0.6});
  expect_probability_zero(lccdf1, 0, 1, {0.3, 0.6}, NEG_INF);
}

TEST(ProbDiscreteCdfsEdgeElements, neg_binomial_2_lcdf_expect_ad) {
  // The infinite count is removed before the beta lcdf call
  auto f = [](const auto& mu, const auto& phi) {
    return stan::math::neg_binomial_2_lcdf(std::vector<int>{3, INT_MAX_N}, mu,
                                           phi);
  };
  Eigen::VectorXd mu(2), phi(2);
  mu << 4.0, 1.5;
  phi << 2.0, 3.0;
  stan::test::expect_ad(f, 4.0, 2.0);
  stan::test::expect_ad(f, mu, 2.0);
  stan::test::expect_ad(f, 4.0, phi);
  stan::test::expect_ad(f, mu, phi);
}
