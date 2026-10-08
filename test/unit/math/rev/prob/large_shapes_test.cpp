#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

// Accuracy at large shape parameters of the distributions whose value or
// partials contain a difference f(x + k) - f(x), f in {lgamma, lbeta,
// digamma}, of a large shape x and a count or offset k. The plain
// differences keep no correct digits once x
// is near 1e15. beta_binomial_lpmf has its own tests in
// test/unit/math/{prim,rev}/prob/beta_binomial*_test.cpp.
//
// References: mpmath at 90 digits (value from mp.loggamma, partials from
// mp.digamma), checked against mp.diff at 140 digits; the cdfs as sums of
// the pmf. The arguments are written in hex so that they are exact. The
// partials of positive parameters are compared as x * d/dx, the gradient in
// log(x) that a sampler sees; the others are compared as they are.
//
// For each function the first rows are points where develop was wrong by
// more than 100 times the tolerance, and the last row is a point at a small
// shape where develop was correct. The cdfs of beta_binomial and
// beta_neg_binomial are limited by the tolerance of their 3F2 series, so
// their tolerances are larger.

namespace large_shapes_test_internal {

struct TestCase {
  const char* tag;
  std::vector<double> args;
  double value;
  std::vector<double> grads;
};

// (value, gradient) tolerance factors, applied as tol * max(1, |reference|)
std::pair<double, double> tolerances(const std::string& tag) {
  if (tag == "BBC" || tag == "BBLC" || tag == "BBLCC" || tag == "BNBC"
      || tag == "BNBLCC") {
    return {1e-6, 1e-4};
  }
  return {1e-12, 1e-11};
}

// for each gradient, the index of the argument it is taken with respect
// to if that argument is a positive parameter (compared on the log scale),
// or -1
std::vector<int> log_scale_args(const std::string& tag) {
  if (tag == "NB") {
    return {1, 2};
  } else if (tag == "NB2") {
    return {-1, 2};
  } else if (tag == "NB2L") {
    return {-1, 2};
  } else if (tag == "GLM1") {
    return {-1, -1, 4};
  } else if (tag == "GLM5") {
    return {-1, -1, 12};
  } else if (tag == "BNB" || tag == "BNBC" || tag == "BNBLCC") {
    return {1, 2, 3};
  } else if (tag == "DM") {
    return {3, 4, 5};
  } else if (tag == "LKJ") {
    return {3};
  } else if (tag == "LKJC") {
    return {5};
  } else if (tag == "BBC" || tag == "BBLC" || tag == "BBLCC") {
    return {2, 3};
  } else if (tag == "ST") {
    return {-1, 1, -1, 3};
  } else if (tag == "MST" || tag == "MSTF") {
    return {-1, -1, 2};
  } else if (tag == "LBETA") {
    return {0, 1};
  }
  return {1};  // YS, YSC, YSLC, YSLCC
}

void eval(const std::string& tag, const std::vector<double>& a, double& value,
          std::vector<double>& grads) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;
  var lp;
  std::vector<var> wrt;
  auto to_int = [](double x) { return static_cast<int>(x); };
  if (tag == "NB" || tag == "NB2" || tag == "NB2L") {
    var p1 = a[1], p2 = a[2];
    if (tag == "NB") {
      lp = stan::math::neg_binomial_lpmf(to_int(a[0]), p1, p2);
    } else if (tag == "NB2") {
      lp = stan::math::neg_binomial_2_lpmf(to_int(a[0]), p1, p2);
    } else {
      lp = stan::math::neg_binomial_2_log_lpmf(to_int(a[0]), p1, p2);
    }
    wrt = {p1, p2};
  } else if (tag == "GLM1" || tag == "GLM5") {
    const int n_obs = (tag == "GLM1") ? 1 : 5;
    std::vector<int> y(n_obs);
    Eigen::MatrixXd x(n_obs, 1);
    for (int i = 0; i < n_obs; ++i) {
      y[i] = to_int(a[i]);
      x(i, 0) = a[n_obs + i];
    }
    var alpha = a[2 * n_obs];
    Matrix<var, Dynamic, 1> beta(1);
    beta(0) = a[2 * n_obs + 1];
    var phi = a[2 * n_obs + 2];
    lp = stan::math::neg_binomial_2_log_glm_lpmf(y, x, alpha, beta, phi);
    wrt = {alpha, beta(0), phi};
  } else if (tag == "BNB" || tag == "BNBC" || tag == "BNBLCC") {
    var r = a[1], alpha = a[2], beta = a[3];
    const int n = to_int(a[0]);
    if (tag == "BNB") {
      lp = stan::math::beta_neg_binomial_lpmf(n, r, alpha, beta);
    } else if (tag == "BNBC") {
      lp = stan::math::beta_neg_binomial_cdf(n, r, alpha, beta);
    } else {
      lp = stan::math::beta_neg_binomial_lccdf(n, r, alpha, beta);
    }
    wrt = {r, alpha, beta};
  } else if (tag == "DM") {
    std::vector<int> ns{to_int(a[0]), to_int(a[1]), to_int(a[2])};
    Matrix<var, Dynamic, 1> alpha(3);
    alpha << a[3], a[4], a[5];
    lp = stan::math::dirichlet_multinomial_lpmf(ns, alpha);
    wrt = {alpha(0), alpha(1), alpha(2)};
  } else if (tag == "LKJ") {
    Eigen::MatrixXd omega(3, 3);
    omega << 1.0, a[0], a[1], a[0], 1.0, a[2], a[1], a[2], 1.0;
    var eta = a[3];
    lp = stan::math::lkj_corr_lpdf(omega, eta);
    wrt = {eta};
  } else if (tag == "LKJC") {
    Eigen::MatrixXd L = Eigen::MatrixXd::Zero(3, 3);
    L(0, 0) = 1.0;
    L(1, 0) = a[0];
    L(1, 1) = a[1];
    L(2, 0) = a[2];
    L(2, 1) = a[3];
    L(2, 2) = a[4];
    var eta = a[5];
    lp = stan::math::lkj_corr_cholesky_lpdf(L, eta);
    wrt = {eta};
  } else if (tag == "BBC" || tag == "BBLC" || tag == "BBLCC") {
    var alpha = a[2], beta = a[3];
    const int n = to_int(a[0]), N = to_int(a[1]);
    if (tag == "BBC") {
      lp = stan::math::beta_binomial_cdf(n, N, alpha, beta);
    } else if (tag == "BBLC") {
      lp = stan::math::beta_binomial_lcdf(n, N, alpha, beta);
    } else {
      lp = stan::math::beta_binomial_lccdf(n, N, alpha, beta);
    }
    wrt = {alpha, beta};
  } else if (tag == "ST") {
    var y = a[0], nu = a[1], mu = a[2], sigma = a[3];
    lp = stan::math::student_t_lpdf(y, nu, mu, sigma);
    wrt = {y, nu, mu, sigma};
  } else if (tag == "MST" || tag == "MSTF") {
    // dimension 2, mu = 0, identity scale: y1, y2, nu
    Matrix<var, Dynamic, 1> y(2);
    y << a[0], a[1];
    Eigen::VectorXd mu = Eigen::VectorXd::Zero(2);
    Eigen::MatrixXd scale = Eigen::MatrixXd::Identity(2, 2);
    var nu = a[2];
    if (tag == "MST") {
      lp = stan::math::multi_student_t_cholesky_lpdf(y, nu, mu, scale);
    } else {
      lp = stan::math::multi_student_t_lpdf(y, nu, mu, scale);
    }
    wrt = {y(0), y(1), nu};
  } else if (tag == "LBETA") {
    var d = a[0], x = a[1];
    lp = stan::math::lbeta(d, x);
    wrt = {d, x};
  } else if (tag == "YS") {
    var alpha = a[1];
    lp = stan::math::yule_simon_lpmf(to_int(a[0]), alpha);
    wrt = {alpha};
  } else {
    var alpha = a[1];
    const int n = to_int(a[0]);
    if (tag == "YSC") {
      lp = stan::math::yule_simon_cdf(n, alpha);
    } else if (tag == "YSLC") {
      lp = stan::math::yule_simon_lcdf(n, alpha);
    } else {
      lp = stan::math::yule_simon_lccdf(n, alpha);
    }
    wrt = {alpha};
  }
  lp.grad();
  value = lp.val();
  grads.clear();
  for (const var& v : wrt) {
    grads.push_back(v.adj());
  }
  stan::math::recover_memory();
}

// clang-format off
const std::vector<TestCase> test_cases = {
    // neg_binomial_lpmf(n | alpha, beta): n, alpha, beta
    {"NB", {0x1.8000000000000p+1, 0x1.c6bf526340000p+49, 0x1.2f2a36ecd5555p+48}, -1.4959226032237274, {1.3124999999999966e-30, 5.6249999999999838e-31}},
    {"NB", {0x1.0000000000000p+0, 0x1.c6bf526340000p+49, 0x1.2f2a36ecd5555p+48}, -1.9013877113318889, {-1.9999999999999957e-15, 5.9999999999999829e-15}},
    {"NB", {0x1.9000000000000p+5, 0x1.c6bf526340000p+49, 0x1.2309ce5400000p+44}, -2.8766166803657541, {2.4999999999998758e-29, 0.0}},
    {"NB", {0x0.0p+0, 0x1.9000000000000p+6, 0x1.0aaaaaaaaaaabp+5}, -2.9558802241544401, {-2.9558802241544401e-2, 8.7378640776699017e-2}},
    // neg_binomial_2_lpmf(n | mu, phi): n, mu, phi
    {"NB2", {0x1.c800000000000p+5, 0x1.8000000000000p+1, 0x1.2a05f20000000p+33}, -1.1677494780996510e+2, {1.7999999994600000e+1, -1.4294999940379000e-17}},
    {"NB2", {0x1.4000000000000p+2, 0x1.9000000000000p+5, 0x1.2a05f20000000p+33}, -3.5227376614641316e+1, {-8.9999999550000002e-1, -1.0099999929136667e-17}},
    {"NB2", {0x1.0000000000000p+0, 0x1.8000000000000p+1, 0x1.2a05f20000000p+33}, -1.9013877111818903, {-6.6666666646666667e-1, -1.4999999991000000e-20}},
    {"NB2", {0x0.0p+0, 0x1.8000000000000p+1, 0x1.9000000000000p+6}, -2.9558802241544403, {-9.7087378640776699e-1, -4.3258864931139302e-4}},
    // neg_binomial_2_log_lpmf(n | eta, phi): n, eta, phi
    {"NB2L", {0x1.8000000000000p+1, 0x1.193ea7aad030bp+0, 0x1.c6bf526340000p+49}, -1.4959226032237274, {-2.7213891705004510e-16, 1.4999999999999960e-30}},
    {"NB2L", {0x1.4000000000000p+2, 0x1.f4bd2b7ac1bafp+1, 0x1.c6bf526340000p+49}, -3.5227376715640301e+1, {-4.4999999999997745e+1, -1.0099999999999289e-27}},
    {"NB2L", {0x1.0000000000000p+1, -0x1.6d3c324e13f50p-2, 0x1.37807ed5e8000p+50}, -2.1064970684374103, {1.2999999999999994, 8.2582982577654707e-32}},
    {"NB2L", {0x0.0p+0, 0x1.193ea7aad030bp+0, 0x1.9000000000000p+6}, -2.9558802241544405, {-2.9126213592233012, -4.3258864931139310e-4}},
    // neg_binomial_2_log_glm_lpmf, one observation: y, x, alpha, beta, phi
    {"GLM1", {0x0.0p+0, 0x1.0000000000000p-1, 0x1.0000000000000p+0, 0x1.8000000000000p-1, 0x1.c6bf526340000p+49}, -3.9550767229205693, {-3.9550767229205615, -1.9775383614602807, -7.8213159420940446e-30}},
    {"GLM1", {0x1.4000000000000p+2, 0x1.0000000000000p-1, 0x1.0000000000000p+0, 0x1.8000000000000p-1, 0x1.c6bf526340000p+49}, -1.8675684657026251, {1.0449232770794187, 5.2246163853970937e-1, 1.9540676725087929e-30}},
    {"GLM1", {0x0.0p+0, 0x1.0000000000000p-1, -0x1.999999999999ap-3, 0x1.8000000000000p-1, 0x1.37807ed5e8000p+50}, -1.1912462166123576, {-1.1912462166123571, -5.9562310830617854e-1, -3.7803493755481261e-31}},
    {"GLM1", {0x0.0p+0, 0x1.0000000000000p-1, 0x1.0000000000000p+0, 0x1.8000000000000p-1, 0x1.9000000000000p+6}, -3.8788665246722319, {-3.8046018026251338, -1.9023009013125669, -7.4264722047098149e-4}},
    // neg_binomial_2_log_glm_lpmf, five observations: y1..y5, x1..x5, alpha, beta, phi
    {"GLM5", {0x1.0000000000000p+1, 0x1.0000000000000p+2, 0x1.8000000000000p+1, 0x1.4000000000000p+2, 0x1.8000000000000p+2, 0x1.0000000000000p-1, -0x1.0000000000000p-2, 0x1.4000000000000p+0, 0x0.0p+0, 0x1.0000000000000p+1, 0x1.0000000000000p+0, 0x1.8000000000000p-1, 0x1.c6bf526340000p+49}, -1.3267966555421410e+1, {-8.0507631204932393, -1.8705862362560038e+1, -2.2918189242185528e-29}},
    {"GLM5", {0x0.0p+0, 0x1.0000000000000p+0, 0x1.4000000000000p+2, 0x1.8000000000000p+1, 0x1.c800000000000p+5, 0x1.0000000000000p-1, -0x1.0000000000000p-2, 0x1.4000000000000p+0, 0x0.0p+0, 0x1.0000000000000p+1, 0x1.0000000000000p+0, 0x1.8000000000000p-1, 0x1.c6bf526340000p+49}, -5.5025862739499810e+1, {3.7949236879506146e+1, 8.5544137637438705e+1, -9.8183556706826074e-28}},
    {"GLM5", {0x1.0000000000000p+1, 0x1.0000000000000p+2, 0x1.8000000000000p+1, 0x1.4000000000000p+2, 0x1.8000000000000p+2, 0x1.0000000000000p-1, -0x1.0000000000000p-2, 0x1.4000000000000p+0, 0x0.0p+0, 0x1.0000000000000p+1, 0x1.0000000000000p+0, 0x1.8000000000000p-1, 0x1.d1a94a2000000p+39}, -1.3267966555398514e+1, {-8.0507631203930684, -1.8705862362370543e+1, -2.2918189241713975e-23}},
    {"GLM5", {0x0.0p+0, 0x1.0000000000000p+0, 0x1.4000000000000p+2, 0x1.8000000000000p+1, 0x1.c800000000000p+5, 0x1.0000000000000p-1, -0x1.0000000000000p-2, 0x1.4000000000000p+0, 0x0.0p+0, 0x1.0000000000000p+1, 0x1.0000000000000p+0, 0x1.8000000000000p-1, 0x1.9000000000000p+6}, -4.7240731087544224e+1, {3.3378922648644250e+1, 7.6036039699167587e+1, -6.2120571346832957e-2}},
    // beta_neg_binomial_lpmf(n | r, alpha, beta): n, r, alpha, beta
    {"BNB", {0x1.4000000000000p+2, 0x1.4000000000000p+2, 0x1.c6bf526340000p+49, 0x1.c6bf526340000p+50}, -2.6841050769298924, {-3.5297736803318772e-1, 1.6666666666666617e-15, -8.3333333333333083e-16}},
    {"BNB", {0x1.4000000000000p+3, 0x1.4000000000000p+2, 0x1.c6bf526340000p+49, 0x1.c6bf526340000p+50}, -2.6389577451069742, {6.9616704560883204e-2, 1.6666666666666591e-30, 4.1666666666666470e-31}},
    {"BNB", {0x1.0000000000000p+0, 0x1.4000000000000p+2, 0x1.c6bf526340000p+49, 0x1.c6bf526340000p+50}, -4.2890886390146075, {-8.9861228866810702e-1, 2.9999999999999917e-15, -1.4999999999999983e-15}},
    {"BNB", {0x0.0p+0, 0x1.9000000000000p+6, 0x1.9000000000000p+6, 0x1.9000000000000p+6}, -5.2468933188272313e+1, {-4.0629959884472565e-1, 2.8935383163709856e-1, -4.0629959884472565e-1}},
    // dirichlet_multinomial_lpmf(ns | alpha): n1, n2, n3, alpha1, alpha2, alpha3
    {"DM", {0x0.0p+0, 0x1.0000000000000p+0, 0x0.0p+0, 0x1.6bcc41e900000p+47, 0x1.10d9316ec0000p+48, 0x1.c6bf526340000p+48}, -1.2039728043259360, {-1.0000000000000000e-15, 2.3333333333333333e-15, -1.0000000000000000e-15}},
    {"DM", {0x1.0000000000000p+0, 0x1.0000000000000p+1, 0x1.0000000000000p+1, 0x1.6bcc41e900000p+47, 0x1.10d9316ec0000p+48, 0x1.c6bf526340000p+48}, -2.0024805005437123, {9.9999999999999700e-30, 1.6666666666666656e-15, -9.9999999999999400e-16}},
    {"DM", {0x1.4000000000000p+2, 0x0.0p+0, 0x0.0p+0, 0x1.6bcc41e900000p+47, 0x1.10d9316ec0000p+48, 0x1.c6bf526340000p+48}, -8.0471895621704619, {1.9999999999999760e-14, -4.9999999999999900e-15, -4.9999999999999900e-15}},
    {"DM", {0x0.0p+0, 0x1.0000000000000p+0, 0x0.0p+0, 0x1.4000000000000p+4, 0x1.e000000000000p+4, 0x1.9000000000000p+5}, -1.2039728043259360, {-1.0000000000000000e-2, 2.3333333333333333e-2, -1.0000000000000000e-2}},
    // lkj_corr_lpdf(Omega | eta): Omega(1, 0), Omega(2, 0), Omega(2, 1), eta
    {"LKJ", {0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x1.c6bf526340000p+49}, 5.0091069763591928e+1, {1.4999999999999999e-15}},
    {"LKJ", {0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x1.772aa3f848000p+51}, 5.1881953466300579e+1, {4.5454545454545453e-16}},
    {"LKJ", {0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x1.802ba9f400000p+41}, 4.1520320547827412e+1, {4.5454545454544307e-13}},
    {"LKJ", {0x1.3333333333333p-2, -0x1.999999999999ap-3, 0x1.999999999999ap-4, 0x1.9000000000000p+6}, -1.1130679230833304e+1, {-1.4988714303399179e-1}},
    // lkj_corr_cholesky_lpdf(L | eta): L(1, 0), L(1, 1), L(2, 0), L(2, 1), L(2, 2), eta
    {"LKJC", {0x0.0p+0, 0x1.0000000000000p+0, 0x0.0p+0, 0x0.0p+0, 0x1.0000000000000p+0, 0x1.c6bf526340000p+49}, 5.0091069763591928e+1, {1.4999999999999999e-15}},
    {"LKJC", {0x0.0p+0, 0x1.0000000000000p+0, 0x0.0p+0, 0x0.0p+0, 0x1.0000000000000p+0, 0x1.772aa3f848000p+51}, 5.1881953466300579e+1, {4.5454545454545453e-16}},
    {"LKJC", {0x0.0p+0, 0x1.0000000000000p+0, 0x0.0p+0, 0x0.0p+0, 0x1.0000000000000p+0, 0x1.802ba9f400000p+41}, 4.1520320547827412e+1, {4.5454545454544307e-13}},
    {"LKJC", {0x1.3333333333333p-2, 0x1.e86ab810ea912p-1, -0x1.999999999999ap-3, 0x1.57808173fc2ddp-3, 0x1.ee40264218109p-1, 0x1.9000000000000p+6}, -1.1177834570568931e+1, {-1.4988714303399185e-1}},
    // yule_simon_{cdf,lcdf,lccdf}(n | alpha): n, alpha
    {"YSC", {0x1.0000000000000p+0, 0x1.9000000000000p+6}, 9.9009900990099010e-1, {9.8029604940692089e-5}},
    {"YSLC", {0x1.0000000000000p+0, 0x1.9000000000000p+6}, -9.9503308531680828e-3, {9.9009900990099010e-5}},
    {"YSLCC", {0x1.0000000000000p+1, 0x1.c6bf526340000p+49}, -6.8384405609261428e+1, {-1.9999999999999970e-15}},
    {"YSLCC", {0x1.0000000000000p+0, 0x1.c6bf526340000p+49}, -3.4538776394910686e+1, {-9.9999999999999900e-16}},
    {"YSLCC", {0x1.4000000000000p+2, 0x1.c6bf526340000p+49}, -1.6790639023177140e+2, {-4.9999999999999850e-15}},
    {"YSLCC", {0x1.0000000000000p+0, 0x1.9000000000000p+6}, -4.6151205168412595, {-9.9009900990099010e-3}},
    // beta_binomial_{cdf,lcdf,lccdf}(n | N, alpha, beta): n, N, alpha, beta
    {"BBC", {0x1.d000000000000p+4, 0x1.d400000000000p+6, 0x1.c6bf526340000p+49, 0x1.550f7dca70000p+51}, 5.2836214512816302e-1, {-1.8707495486824231e-15, 6.2358318289414104e-16}},
    {"BBC", {0x1.4000000000000p+2, 0x1.d400000000000p+6, 0x1.c6bf526340000p+49, 0x1.550f7dca70000p+51}, 1.9080827003311692e-9, {-4.6543958285535218e-23, 1.5514652761844824e-23}},
    {"BBC", {0x1.0000000000000p+0, 0x1.d400000000000p+6, 0x1.c6bf526340000p+49, 0x1.550f7dca70000p+51}, 9.6433473051139671e-14, {-2.7266564505209334e-27, 9.0888548350696083e-28}},
    {"BBC", {0x0.0p+0, 0x1.d400000000000p+6, 0x1.9000000000000p+6, 0x1.9000000000000p+6}, 2.0588830493356332e-26, {-9.5019176900109660e-27, 6.5044482281184328e-27}},
    {"BBLC", {0x0.0p+0, 0x1.d400000000000p+6, 0x1.3880000000000p+13, 0x1.3880000000000p+13}, -8.0760883216610174e+1, {-5.8331005942734647e-3, 5.7995618892517663e-3}},
    {"BBLC", {0x1.4000000000000p+2, 0x1.d400000000000p+6, 0x1.3880000000000p+13, 0x1.3880000000000p+13}, -6.1834457741779739e+1, {-5.3377886925755501e-3, 5.3097367250871483e-3}},
    {"BBLC", {0x1.4000000000000p+2, 0x1.d400000000000p+6, 0x1.e848000000000p+19, 0x1.e848000000000p+19}, -6.2113769736395483e+1, {-5.3543714610347851e-5, 5.3540877034666121e-5}},
    {"BBLC", {0x1.c800000000000p+5, 0x1.d400000000000p+6, 0x1.9000000000000p+6, 0x1.9000000000000p+6}, -8.1708867859107758e-1, {-3.8676835548171366e-2, 3.8433660832149159e-2}},
    {"BBLCC", {0x1.d000000000000p+4, 0x1.d400000000000p+6, 0x1.c6bf526340000p+49, 0x1.550f7dca70000p+51}, -7.5154384451605549e-1, {3.9664957538041156e-15, -1.3221652512680385e-15}},
    {"BBLCC", {0x1.4000000000000p+2, 0x1.d400000000000p+6, 0x1.c6bf526340000p+49, 0x1.550f7dca70000p+51}, -1.9080827021515590e-9, {4.6543958374344940e-23, -1.5514652791448065e-23}},
    {"BBLCC", {0x1.0000000000000p+0, 0x1.d400000000000p+6, 0x1.c6bf526340000p+49, 0x1.550f7dca70000p+51}, -9.6433473051144321e-14, {2.7266564505211963e-27, -9.0888548350704848e-28}},
    {"BBLCC", {0x0.0p+0, 0x1.d400000000000p+6, 0x1.9000000000000p+6, 0x1.9000000000000p+6}, -2.0588830493356332e-26, {9.5019176900109660e-27, -6.5044482281184328e-27}},
    // beta_neg_binomial_{cdf,lccdf}(n | r, alpha, beta): n, r, alpha, beta
    {"BNBC", {0x0.0p+0, 0x1.4000000000000p+2, 0x1.c6bf526340000p+49, 0x1.c6bf526340000p+50}, 4.1152263374485871e-3, {-4.5210382249716626e-3, 1.3717421124828587e-17, -6.8587105624143073e-18}},
    {"BNBC", {0x1.4000000000000p+2, 0x1.4000000000000p+2, 0x1.c6bf526340000p+49, 0x1.c6bf526340000p+50}, 2.1312808006909550e-1, {-1.1257829350011265e-1, 4.5521516029060565e-16, -2.2760758014530300e-16}},
    {"BNBC", {0x1.0000000000000p+0, 0x1.4000000000000p+2, 0x1.c6bf526340000p+49, 0x1.c6bf526340000p+50}, 1.7832647462277188e-2, {-1.6847681416578131e-2, 5.4869684499314275e-17, -2.7434842249657186e-17}},
    {"BNBC", {0x0.0p+0, 0x1.4000000000000p+3, 0x1.4000000000000p+3, 0x1.4000000000000p+3}, 4.6119797244235025e-3, {-1.9089636233770314e-3, 1.4059955145634730e-3, -1.9089636233770314e-3}},
    {"BNBLCC", {0x0.0p+0, 0x1.4000000000000p+2, 0x1.c6bf526340000p+49, 0x1.c6bf526340000p+50}, -4.1237171838620870e-3, {4.5397202011079093e-3, -1.3774104683195648e-17, 6.8870523415978377e-18}},
    {"BNBLCC", {0x1.4000000000000p+2, 0x1.4000000000000p+2, 0x1.c6bf526340000p+49, 0x1.c6bf526340000p+50}, -2.3968978849662940e-1, {1.4307067090410113e-1, -5.7851239669421455e-16, 2.8925619834710749e-16}},
    {"BNBLCC", {0x1.4000000000000p+3, 0x1.4000000000000p+2, 0x1.c6bf526340000p+49, 0x1.c6bf526340000p+50}, -9.0618005875034648e-1, {3.6151768966033925e-1, -1.7679265277287130e-15, 8.8396326386435604e-16}},
    {"BNBLCC", {0x0.0p+0, 0x1.4000000000000p+3, 0x1.4000000000000p+3, 0x1.4000000000000p+3}, -4.6226477159237375e-3, {1.9178085173744892e-3, -1.4125099819608222e-3, 1.9178085173744892e-3}},
    // student_t_lpdf(y | nu, mu, sigma): y, nu, mu, sigma
    {"ST", {0x1.0000000000000p+0, 0x1.6bcc41e900000p+46, 0x0.0p+0, 0x1.0000000000000p+0}, -1.4189385332046777, {-1.0000000000000000, 4.9999999999999833e-29, 1.0000000000000000, 2.6727647100921956e-51}},
    {"ST", {-0x1.cd1e504efb30cp+1, 0x1.7cc4b890abebfp+45, 0x1.64703afcf3380p-4, 0x1.c43477d4ae376p+3}, -3.6014210338578538, {1.8475570774651872e-2, 1.0330563009181428e-28, -1.8475570774651872e-2, -6.5940664388226970e-2}},
    {"ST", {0x1.8000000000000p+1, 0x1.6bcc41e900000p+46, 0x0.0p+0, 0x1.0000000000000p+1}, -2.7370857137646191, {-7.4999999999999062e-1, 1.0937500000001266e-29, 7.4999999999999062e-1, 6.2499999999998594e-1}},
    {"ST", {-0x1.c8dd60359f470p-1, 0x1.da533b967d6d3p-4, -0x1.572a29ad6231cp-1, 0x1.9d9944276e11dp-3}, -1.6062891011707614, {4.5853946438927804, 6.8870497923078714, -4.5853946438927804, 9.0518700772802596e-2}},
    // multi_student_t_cholesky_lpdf (MST) and multi_student_t_lpdf (MSTF),
    // dimension 2, mu = 0, identity scale: y1, y2, nu
    {"MST", {0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.d1a94a2000000p+39}, -2.8378770664103455, {-1.0000000000000000, -1.0000000000000000, 9.9999999999866667e-25}},
    {"MST", {0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.74876e8000000p+36}, -2.8378770664193455, {-1.0000000000000000, -1.0000000000000000, 9.9999999998666667e-23}},
    {"MST", {0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.dcd6500000000p+29}, -2.8378770674093455, {-1.0000000000000000, -1.0000000000000000, 9.9999999866666667e-19}},
    {"MST", {-0x1.03f559f1746d0p-1, 0x1.22539ac62abe0p-3, 0x1.681d9f093b62ap-6}, -4.4798146609821203, {3.4235928669135459, -9.5588370652474091e-1, 4.1318417102814693e+1}},
    {"MSTF", {0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.d1a94a2000000p+39}, -2.8378770664103455, {-1.0000000000000000, -1.0000000000000000, 9.9999999999866667e-25}},
    {"MSTF", {0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.74876e8000000p+36}, -2.8378770664193455, {-1.0000000000000000, -1.0000000000000000, 9.9999999998666667e-23}},
    {"MSTF", {0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.dcd6500000000p+29}, -2.8378770674093455, {-1.0000000000000000, -1.0000000000000000, 9.9999999866666667e-19}},
    {"MSTF", {-0x1.03f559f1746d0p-1, 0x1.22539ac62abe0p-3, 0x1.681d9f093b62ap-6}, -4.4798146609821203, {3.4235928669135459, -9.5588370652474091e-1, 4.1318417102814693e+1}},
    // yule_simon_lpmf(n | alpha): n, alpha
    {"YS", {0x1.0000000000000p+0, 0x1.d1a94a2000000p+39}, -9.9999999999950000e-13, {9.9999999999900000e-25}},
    {"YS", {0x1.0000000000000p+0, 0x1.74876e8000000p+36}, -9.9999999999500000e-12, {9.9999999999000000e-23}},
    {"YS", {0x1.0000000000000p+0, 0x1.2a05f20000000p+33}, -9.9999999995000000e-11, {9.9999999990000000e-21}},
    {"YS", {0x1.0000000000000p+0, 0x1.cc271540e654dp-2}, -1.1710409729900522, {1.5353924179257093}},
    // lbeta(d, x) with var arguments: the derivative digamma(x) - digamma(d + x)
    {"LBETA", {0x1.8000000000000p+0, 0x1.6bcc41e900000p+46}, -4.8475069190510208e+1, {-3.2199701327938073e+1, -1.4999999999999962e-14}},
    {"LBETA", {0x1.0000000000000p+0, 0x1.6bcc41e900000p+46}, -3.2236191301916640e+1, {-3.2813406966818177e+1, -1.0000000000000000e-14}},
    {"LBETA", {0x1.0000000000000p-1, 0x1.6bcc41e900000p+46}, -1.5545730708033618e+1, {-3.4199701327938063e+1, -5.0000000000000125e-15}},
    {"LBETA", {0x1.ba1a7a7909c76p-7, 0x1.f4f68d902030cp-6}, 4.6705195755596318, {-5.1474613897767773e+1, -1.0033990131067132e+1}},
};
// clang-format on
}  // namespace large_shapes_test_internal

TEST(ProbDistributions, large_shapes_value_and_log_scale_gradients) {
  using large_shapes_test_internal::eval;
  using large_shapes_test_internal::log_scale_args;
  using large_shapes_test_internal::test_cases;
  using large_shapes_test_internal::tolerances;
  for (const auto& t : test_cases) {
    const std::string tag(t.tag);
    double value;
    std::vector<double> grads;
    eval(tag, t.args, value, grads);
    const auto tol = tolerances(tag);
    EXPECT_NEAR(value, t.value, tol.first * std::max(1.0, std::fabs(t.value)))
        << tag << " value, last argument " << t.args.back();
    const std::vector<int> log_args = log_scale_args(tag);
    ASSERT_EQ(grads.size(), t.grads.size()) << tag;
    for (size_t j = 0; j < grads.size(); ++j) {
      const double scale = log_args[j] >= 0 ? t.args[log_args[j]] : 1.0;
      const double got = scale * grads[j];
      const double expected = scale * t.grads[j];
      EXPECT_NEAR(got, expected,
                  tol.second * std::max(1.0, std::fabs(expected)))
          << tag << " gradient " << j << ", last argument " << t.args.back();
    }
  }
}

TEST(ProbDistributions, lkj_corr_eta_gradient_at_one) {
  // At eta == 1.0 exactly, develop returned the constant without its
  // derivative, so d/deta lost sum_k psi(eta + (K - 1) / 2) - psi(eta +
  // (K - 1 - k) / 2). eta = exp(0) = 1 is the value at the default
  // initialization of an unconstrained parameter. Reference from mpmath.
  Eigen::MatrixXd omega(3, 3);
  omega << 1.0, 0.3, -0.2, 0.3, 1.0, 0.1, -0.2, 0.1, 1.0;
  stan::math::var eta = 1.0;
  stan::math::var lp = stan::math::lkj_corr_lpdf(omega, eta);
  lp.grad();
  EXPECT_NEAR(lp.val(), -1.596312591138855, 1e-14);
  EXPECT_NEAR(eta.adj(), 1.2214197179296566, 1e-14);
  stan::math::recover_memory();
}

TEST(ProbDistributions, large_shapes_propto_dropped_terms) {
  // full - propto must be the sum of the terms that propto drops, computed
  // here from their formulas. Some of these terms are inside a function
  // call: lgamma(n) in lbeta(n, alpha + 1), lgamma(y + 1) in lchoose.
  using stan::math::var;
  using std::lgamma;
  using std::log;

  // student_t_lpdf, parameter nu (both forms of the nu terms: nu < 40 and
  // nu >= 40), data y, mu, and sigma = 1: drops log(sqrt(pi)) per term
  {
    Eigen::VectorXd y(3);
    y << 0.3, -1.2, 2.0;
    Eigen::VectorXd nu_d(3);
    nu_d << 3.5, 55.0, 500.0;
    Eigen::Matrix<var, Eigen::Dynamic, 1> nu = nu_d.cast<var>();
    var full = stan::math::student_t_lpdf<false>(y, nu, 0.1, 1.0);
    var prop = stan::math::student_t_lpdf<true>(y, nu, 0.1, 1.0);
    EXPECT_NEAR(full.val() - prop.val(), -3 * stan::math::LOG_SQRT_PI, 1e-13);
    stan::math::recover_memory();
  }
  // yule_simon_lpmf, parameter alpha: drops lgamma(n)
  {
    std::vector<int> n{1, 3, 7};
    Eigen::VectorXd alpha_d(3);
    alpha_d << 2.5, 0.7, 30.0;
    Eigen::Matrix<var, Eigen::Dynamic, 1> alpha = alpha_d.cast<var>();
    var full = stan::math::yule_simon_lpmf<false>(n, alpha);
    var prop = stan::math::yule_simon_lpmf<true>(n, alpha);
    EXPECT_NEAR(full.val() - prop.val(),
                lgamma(1.0) + lgamma(3.0) + lgamma(7.0), 1e-13);
    stan::math::recover_memory();
  }
  // beta_neg_binomial_lpmf, parameters r, alpha, beta: drops -lgamma(n + 1)
  {
    std::vector<int> n{3, 5};
    var r = 2.5;
    var a = 3.5;
    var b = 1.5;
    var full = stan::math::beta_neg_binomial_lpmf<false>(n, r, a, b);
    var prop = stan::math::beta_neg_binomial_lpmf<true>(n, r, a, b);
    EXPECT_NEAR(full.val() - prop.val(), -(lgamma(4.0) + lgamma(6.0)), 1e-13);
    stan::math::recover_memory();
  }
  // neg_binomial_2_log_glm_lpmf: drops -lgamma(y + 1); y theta if x, alpha
  // and beta are data; phi log(phi) - lgamma(phi) + lgamma(y + phi) if phi
  // is data
  {
    std::vector<int> y{3, 5, 2};
    Eigen::MatrixXd x(3, 2);
    x << 0.5, -0.2, 1.1, 0.3, -0.4, 0.8;
    Eigen::VectorXd beta_d(2);
    beta_d << 0.4, -0.3;
    const double alpha = 0.3;
    const double phi_d = 3.5;
    Eigen::VectorXd theta = (x * beta_d).array() + alpha;
    double lgamma_y1 = 0;
    double y_theta = 0;
    double phi_terms = 0;
    for (size_t i = 0; i < y.size(); ++i) {
      lgamma_y1 += lgamma(y[i] + 1.0);
      y_theta += y[i] * theta(i);
      phi_terms += phi_d * log(phi_d) - lgamma(phi_d) + lgamma(y[i] + phi_d);
    }
    Eigen::Matrix<var, Eigen::Dynamic, 1> beta = beta_d.cast<var>();
    var phi = phi_d;
    // parameters beta and phi
    var full = stan::math::neg_binomial_2_log_glm_lpmf<false>(y, x, alpha, beta,
                                                              phi);
    var prop
        = stan::math::neg_binomial_2_log_glm_lpmf<true>(y, x, alpha, beta, phi);
    EXPECT_NEAR(full.val() - prop.val(), -lgamma_y1, 1e-13);
    // parameter phi; x, alpha and beta are data
    full = stan::math::neg_binomial_2_log_glm_lpmf<false>(y, x, alpha, beta_d,
                                                          phi);
    prop = stan::math::neg_binomial_2_log_glm_lpmf<true>(y, x, alpha, beta_d,
                                                         phi);
    EXPECT_NEAR(full.val() - prop.val(), -lgamma_y1 + y_theta, 1e-13);
    // parameter beta; phi is data
    full = stan::math::neg_binomial_2_log_glm_lpmf<false>(y, x, alpha, beta,
                                                          phi_d);
    prop = stan::math::neg_binomial_2_log_glm_lpmf<true>(y, x, alpha, beta,
                                                         phi_d);
    EXPECT_NEAR(full.val() - prop.val(), -lgamma_y1 + phi_terms, 1e-13);
    stan::math::recover_memory();
  }
}
