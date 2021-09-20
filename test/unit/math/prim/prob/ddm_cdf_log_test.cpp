#include <stan/math/prim/prob/ddm_lcdf.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbDdm, ddm_lcdf_matches_known_cdf) {
  using stan::math::ddm_lcdf;
  using std::exp;
  using std::vector;

  // Note:
  // 1. Each expected log_prob is calculated with the R package `fddm` using
  //    fddm::pfddm(rt, response, a, v, t0, w, sv, log = TRUE)
  // 2. We must define a large error tolerance because we must check on the log
  //    scale, and since the sum of the log PDFs is very negative (~ -1800),
  //    exponentiating this very negative value would result in rounding to
  //    zero. This value of error tolerance is slightly arbitrary, but it is
  //    more useful than comparing zero to zero due to rounding issues.

  static const double pfddm_output = -972.3893812;
  static const double err_tol = 1.0;
  static const vector<double> rt{0.1, 1, 10.0};
  static const int response = 1;  // "lower" threshold
  static const vector<double> a{0.5, 1.0, 5.0};
  static const vector<double> v{-2.0, 0.0, 2.0};
  static const double t0 = 0.0001;
  static const vector<double> w{0.2, 0.5, 0.8};
  static const vector<double> sv{0.0, 0.5, 1.0, 1.5};

  static const int n_rt = rt.size();
  static const int n_a = a.size();
  static const int n_v = v.size();
  static const int n_w = w.size();
  static const int n_sv = sv.size();
  static const int n_max = n_rt * n_a * n_v * n_w * n_sv;
  vector<double> rt_vec(n_max);
  vector<double> a_vec(n_max);
  vector<double> v_vec(n_max);
  vector<double> t0_vec(n_max, t0);
  vector<double> w_vec(n_max);
  vector<double> sv_vec(n_max);

  for (int i = 0; i < n_rt; i++) {
    for (int j = i; j < n_max; j += n_rt) {
      rt_vec[j] = rt[i];
    }
  }
  for (int i = 0; i < n_a; i++) {
    for (int j = i; j < n_max; j += n_a) {
      a_vec[j] = a[i];
    }
  }
  for (int i = 0; i < n_v; i++) {
    for (int j = i; j < n_max; j += n_v) {
      v_vec[j] = v[i];
    }
  }
  for (int i = 0; i < n_w; i++) {
    for (int j = i; j < n_max; j += n_w) {
      w_vec[j] = w[i];
    }
  }
  for (int i = 0; i < n_sv; i++) {
    for (int j = i; j < n_max; j += n_sv) {
      sv_vec[j] = sv[i];
    }
  }

  EXPECT_NEAR((ddm_lcdf(rt, response, a, v, t0_vec, w, sv_vec)), pfddm_output,
              err_tol);
  EXPECT_NEAR((ddm_lcdf<true>(rt, response, a, v, t0_vec, w, sv_vec)),
              0,  // true makes ddm_lcdf() and wiener_lpdf() evaluate to 0
              err_tol);
  EXPECT_NEAR((ddm_lcdf<false>(rt, response, a, v, t0_vec, w, sv_vec)),
              pfddm_output, err_tol);
  EXPECT_NEAR((ddm_lcdf<vector<double>, int, vector<double>, vector<double>,
                        vector<double>, vector<double>, vector<double> >(
                  rt, response, a, v, t0_vec, w, sv_vec)),
              pfddm_output, err_tol);
  EXPECT_NEAR(
      (ddm_lcdf<true, vector<double>, int, vector<double>, vector<double>,
                vector<double>, vector<double>, vector<double> >(
          rt, response, a, v, t0_vec, w, sv_vec)),
      0,  // true makes ddm_lcdf() and wiener_lpdf() evaluate to 0
      err_tol);
  EXPECT_NEAR(
      (ddm_lcdf<false, vector<double>, int, vector<double>, vector<double>,
                vector<double>, vector<double>, vector<double> >(
          rt, response, a, v, t0_vec, w, sv_vec)),
      pfddm_output, err_tol);
}
