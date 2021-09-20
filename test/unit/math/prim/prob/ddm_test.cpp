#include <stan/math/prim/prob/ddm_lpdf.hpp>
#include <stan/math/prim/prob/wiener_lpdf.hpp>
#include <gtest/gtest.h>
#include <vector>

// ddm_lpdf(rt, response, a, v, t0, w, sv)
// wiener_lpdf(y, alpha, tau, beta, delta)
// rt <- y
// response <- 2 (wiener_lpdf() always uses the "upper" threshold)
// a <- alpha
// v <- delta
// t0 <- tau
// w <- beta
// sv is not included in wiener_lpdf()
// alpha -> a
// tau -> t0
// beta -> w
// delta -> v
// Note: `response` and `sv` are not included in wiener_lpdf()


// Check invalid arguments

// rt
TEST(mathPrimScalProbDdmScal, invalid_rt) {
  using stan::math::ddm_lpdf;
  using stan::math::INFTY;
  EXPECT_THROW(ddm_lpdf(0, 2, 1, -1, 0, 0.5, 0), std::domain_error);
  EXPECT_THROW(ddm_lpdf(-1, 2, 1, -1, 0, 0.5, 0), std::domain_error);
  EXPECT_THROW(ddm_lpdf(INFTY, 2, 1, -1, 0, 0.5, 0),
               std::domain_error);
  EXPECT_THROW(ddm_lpdf(-INFTY, 2, 1, -1, 0, 0.5, 0),
               std::domain_error);
  EXPECT_THROW(ddm_lpdf(NAN, 2, 1, -1, 0, 0.5, 0), std::domain_error);
}
TEST(mathPrimScalProbDdmMat, invalid_rt) {
  using stan::math::ddm_lpdf;
  using stan::math::INFTY;
  std::vector<double> rt{1, 0};
  EXPECT_THROW(ddm_lpdf(rt, 2, 1, -1, 0, 0.5, 0), std::domain_error);
  rt[1] = -1;
  EXPECT_THROW(ddm_lpdf(rt, 2, 1, -1, 0, 0.5, 0), std::domain_error);
  rt[1] = INFTY;
  EXPECT_THROW(ddm_lpdf(rt, 2, 1, -1, 0, 0.5, 0), std::domain_error);
  rt[1] = -INFTY;
  EXPECT_THROW(ddm_lpdf(rt, 2, 1, -1, 0, 0.5, 0), std::domain_error);
  rt[1] = NAN;
  EXPECT_THROW(ddm_lpdf(rt, 2, 1, -1, 0, 0.5, 0), std::domain_error);
}

// response
TEST(mathPrimScalProbDdmScal, invalid_response) {
  using stan::math::ddm_lpdf;
  using stan::math::INFTY;
  EXPECT_THROW(ddm_lpdf(1, 0, 1, -1, 0, 0.5, 0), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, 3, 1, -1, 0, 0.5, 0), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, -1, 1, -1, 0, 0.5, 0), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, INFTY, 1, -1, 0, 0.5, 0), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, -INFTY, 1, -1, 0, 0.5, 0), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, NAN, 1, -1, 0, 0.5, 0), std::domain_error);
}
TEST(mathPrimScalProbDdmMat, invalid_response) {
  using stan::math::ddm_lpdf;
  using stan::math::INFTY;
  std::vector<double> response{2, 0};
  EXPECT_THROW(ddm_lpdf(1, response, 1, -1, 0, 0.5, 0), std::domain_error);
  response[1] = 3;
  EXPECT_THROW(ddm_lpdf(1, response, 1, -1, 0, 0.5, 0), std::domain_error);
  response[1] = -1;
  EXPECT_THROW(ddm_lpdf(1, response, 1, -1, 0, 0.5, 0), std::domain_error);
  response[1] = INFTY;
  EXPECT_THROW(ddm_lpdf(1, response, 1, -1, 0, 0.5, 0), std::domain_error);
  response[1] = -INFTY;
  EXPECT_THROW(ddm_lpdf(1, response, 1, -1, 0, 0.5, 0), std::domain_error);
  response[1] = NAN;
  EXPECT_THROW(ddm_lpdf(1, response, 1, -1, 0, 0.5, 0), std::domain_error);
}

// a
TEST(mathPrimScalProbDdmScal, invalid_a) {
  using stan::math::ddm_lpdf;
  using stan::math::INFTY;
  EXPECT_THROW(ddm_lpdf(1, 2, 0, -1, 0, 0.5, 0), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, 2, -1, -1, 0, 0.5, 0), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, 2, INFTY, -1, 0, 0.5, 0), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, 2, -INFTY, -1, 0, 0.5, 0), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, 2, NAN, -1, 0, 0.5, 0), std::domain_error);
}
TEST(mathPrimScalProbDdmMat, invalid_a) {
  using stan::math::ddm_lpdf;
  using stan::math::INFTY;
  std::vector<double> a{1, 0};
  EXPECT_THROW(ddm_lpdf(1, 2, a, -1, 0, 0.5, 0), std::domain_error);
  a[1] = -1;
  EXPECT_THROW(ddm_lpdf(1, 2, a, -1, 0, 0.5, 0), std::domain_error);
  a[1] = INFTY;
  EXPECT_THROW(ddm_lpdf(1, 2, a, -1, 0, 0.5, 0), std::domain_error);
  a[1] = -INFTY;
  EXPECT_THROW(ddm_lpdf(1, 2, a, -1, 0, 0.5, 0), std::domain_error);
  a[1] = NAN;
  EXPECT_THROW(ddm_lpdf(1, 2, a, -1, 0, 0.5, 0), std::domain_error);
}

// v
TEST(mathPrimScalProbDdmScal, invalid_v) {
  using stan::math::ddm_lpdf;
  using stan::math::INFTY;
  EXPECT_THROW(ddm_lpdf(1, 2, 1, INFTY, 0, 0.5, 0), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -INFTY, 0, 0.5, 0), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, 2, 1, NAN, 0, 0.5, 0), std::domain_error);
}
TEST(mathPrimScalProbDdmMat, invalid_v) {
  using stan::math::ddm_lpdf;
  using stan::math::INFTY;
  std::vector<double> v{1, INFTY};
  EXPECT_THROW(ddm_lpdf(1, 2, 1, v, 0, 0.5, 0), std::domain_error);
  v[1] = -INFTY;
  EXPECT_THROW(ddm_lpdf(1, 2, 1, v, 0, 0.5, 0), std::domain_error);
  v[1] = NAN;
  EXPECT_THROW(ddm_lpdf(1, 2, 1, v, 0, 0.5, 0), std::domain_error);
}

// t0
TEST(mathPrimScalProbDdmScal, invalid_t0) {
  using stan::math::ddm_lpdf;
  using stan::math::INFTY;
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, -1, 0.5, 0), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, INFTY, 0.5, 0), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, -INFTY, 0.5, 0), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, NAN, 0.5, 0), std::domain_error);
}
TEST(mathPrimScalProbDdmMat, invalid_t0) {
  using stan::math::ddm_lpdf;
  using stan::math::INFTY;
  std::vector<double> t0{1, -1};
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, t0, 0.5, 0), std::domain_error);
  t0[1] = INFTY;
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, t0, 0.5, 0), std::domain_error);
  t0[1] = -INFTY;
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, t0, 0.5, 0), std::domain_error);
  t0[1] = NAN;
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, t0, 0.5, 0), std::domain_error);
}

// w
TEST(mathPrimScalProbDdmScal, invalid_w) {
  using stan::math::ddm_lpdf;
  using stan::math::INFTY;
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, 0, -0.1, 0), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, 0, 0, 0), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, 0, 1, 0), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, 0, 1.1, 0), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, 0, INFTY, 0), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, 0, -INFTY, 0), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, 0, NAN, 0), std::domain_error);
}
TEST(mathPrimScalProbDdmMat, invalid_w) {
  using stan::math::ddm_lpdf;
  using stan::math::INFTY;
  std::vector<double> w{1, -0.1};
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, 0, w, 0), std::domain_error);
  w[1] = 0;
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, 0, w, 0), std::domain_error);
  w[1] = 1;
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, 0, w, 0), std::domain_error);
  w[1] = 1.1;
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, 0, w, 0), std::domain_error);
  w[1] = INFTY;
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, 0, w, 0), std::domain_error);
  w[1] = -INFTY;
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, 0, w, 0), std::domain_error);
  w[1] = NAN;
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, 0, w, 0), std::domain_error);
}

// sv
TEST(mathPrimScalProbDdmScal, invalid_sv) {
  using stan::math::ddm_lpdf;
  using stan::math::INFTY;
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, 0, 0.5, -1), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, 0, 0.5, INFTY), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, 0, 0.5, -INFTY), std::domain_error);
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, 0, 0.5, NAN), std::domain_error);
}
TEST(mathPrimScalProbDdmMat, invalid_sv) {
  using stan::math::ddm_lpdf;
  using stan::math::INFTY;
  std::vector<double> sv{1, -1};
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, 0, 0.5, sv), std::domain_error);
  sv[1] = INFTY;
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, 0, 0.5, sv), std::domain_error);
  sv[1] = -INFTY;
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, 0, 0.5, sv), std::domain_error);
  sv[1] = NAN;
  EXPECT_THROW(ddm_lpdf(1, 2, 1, -1, 0, 0.5, sv), std::domain_error);
}

TEST(ProbDdm, ddm_lpdf_matches_wiener_lpdf) {
  using std::vector;
  using std::exp;
  using stan::math::ddm_lpdf;
  using stan::math::wiener_lpdf;
  // Notes:
  // 1. define error tolerance for PDF approximations, use double tolerance to
  // allow for convergence (of the truncated infinite sum) from above and below
  // 2. wiener_lpdf() only uses the "upper" threshold in the DDM, but the
  // "lower" threshold maps v -> -v and w -> 1-w, and this is covered in the
  // parameter values defined below
  double err_tol = 2 * 0.000001;
  vector<double> rt{0.1, 1, 10.0};
  int response = 2; // wiener_lpdf() always uses the "upper" threshold
  vector<double> a{0.5, 1.0, 5.0};
  vector<double> v{-2.0, 0.0, 2.0};
  double t0 = 0.0001; // t0 (i.e., tau) needs to be > 0 for wiener_lpdf()
  vector<double> w{0.2, 0.5, 0.8};
  double sv = 0.0; // sv is not included in wiener_lpdf(), and thus it must be 0
  
  int n_rt = rt.size();
  int n_a = a.size();
  int n_v = v.size();
  int n_w = w.size();
  int n_max = n_rt * n_a * n_v * n_w;
  vector<double> rt_vec(n_max);
  vector<double> a_vec(n_max);
  vector<double> v_vec(n_max);
  vector<double> t0_vec(n_max, t0);
  vector<double> w_vec(n_max);
  
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
  
  // The PDF approximation error tolerance is based on the non-log version of
  // the PDF. To account for this, we compare the exponentiated version of the
  // *_lpdf results.
  EXPECT_NEAR(exp(ddm_lpdf(rt, response, a, v, t0_vec, w, sv)),
              exp(wiener_lpdf(rt_vec, a_vec, t0_vec, w_vec, v_vec)),
              err_tol);
  EXPECT_NEAR(exp(ddm_lpdf<true>(rt, response, a, v, t0_vec, w, sv)),
              exp(wiener_lpdf<true>(rt_vec, a_vec, t0_vec, w_vec, v_vec)),
              err_tol);
  EXPECT_NEAR(exp(ddm_lpdf<false>(rt, response, a, v, t0_vec, w, sv)),
              exp(wiener_lpdf<false>(rt_vec, a_vec, t0_vec, w_vec, v_vec)),
              err_tol);
  EXPECT_NEAR(
    exp(ddm_lpdf<vector<double>, int, vector<double>, vector<double>,
        vector<double>, vector<double>, double>(
            rt, response, a, v, t0_vec, w, sv)),
            exp(wiener_lpdf<vector<double>, vector<double>, vector<double>,
                vector<double>, vector<double> >(
                    rt_vec, a_vec, t0_vec, w_vec, v_vec)),
                    err_tol);
  EXPECT_NEAR(
    exp(ddm_lpdf<true, vector<double>, int, vector<double>, vector<double>,
        vector<double>, vector<double>, double>(
            rt, response, a, v, t0_vec, w, sv)),
            exp(wiener_lpdf<true, vector<double>, vector<double>, vector<double>,
                vector<double>, vector<double> >(
                    rt_vec, a_vec, t0_vec, w_vec, v_vec)),
                    err_tol);
  EXPECT_NEAR(
    exp(ddm_lpdf<false, vector<double>, int, vector<double>, vector<double>,
        vector<double>, vector<double>, double>(
            rt, response, a, v, t0_vec, w, sv)),
            exp(wiener_lpdf<false, vector<double>, vector<double>, vector<double>,
                vector<double>, vector<double> >(
                    rt_vec, a_vec, t0_vec, w_vec, v_vec)),
                    err_tol);
  
  // check with variable drift rate (against results from R package `fddm`)
  // Notes:
  // 1. We must redefine the error tolerance because we must check on the log
  // scale as the sum of the log PDFs is very negative (~ -1800), and
  // exponentiating this very negative value would result in rounding to zero.
  // This value of error tolerance is slightly arbitrary, but it is more useful
  // than comparing zero to zero due to rounding issues.
  err_tol = 1.0;
  static const double dfddm_output = -1800.6359154;
  vector<double> sv_vals{0.0, 0.5, 1.0, 1.5};
  int n_sv = sv_vals.size();
  n_max *= n_sv;
  vector<double> sv_vec(n_max);
  for (int i = 0; i < n_sv; i++) {
    for (int j = i; j < n_max; j += n_sv) {
      sv_vec[j] = sv_vals[i];
    }
  }
  
  EXPECT_NEAR((ddm_lpdf(rt, response, a, v, t0_vec, w, sv_vec)),
              dfddm_output,
              err_tol);
  EXPECT_NEAR((ddm_lpdf<true>(rt, response, a, v, t0_vec, w, sv_vec)),
              0, // true makes ddm_lpdf() and wiener_lpdf() evaluate to 0
              err_tol);
  EXPECT_NEAR((ddm_lpdf<false>(rt, response, a, v, t0_vec, w, sv_vec)),
              dfddm_output,
              err_tol);
  EXPECT_NEAR(
    (ddm_lpdf<vector<double>, int, vector<double>, vector<double>,
     vector<double>, vector<double>, vector<double> >(
         rt, response, a, v, t0_vec, w, sv_vec)),
         dfddm_output,
         err_tol);
  EXPECT_NEAR(
    (ddm_lpdf<true, vector<double>, int, vector<double>, vector<double>,
     vector<double>, vector<double>, vector<double> >(
         rt, response, a, v, t0_vec, w, sv_vec)),
         0, // true makes ddm_lpdf() and wiener_lpdf() evaluate to 0
         err_tol);
  EXPECT_NEAR(
    (ddm_lpdf<false, vector<double>, int, vector<double>, vector<double>,
     vector<double>, vector<double>, vector<double> >(
         rt, response, a, v, t0_vec, w, sv_vec)),
         dfddm_output,
         err_tol);
}
