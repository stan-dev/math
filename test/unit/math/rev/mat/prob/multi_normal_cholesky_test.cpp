#include <stan/math/prim/mat/fun/cholesky_decompose.hpp>
#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/util.hpp>
// For speed comparisons
//#include <chrono>
//#include <test/unit/math/rev/mat/prob/multi_normal_cholesky_old.hpp>
//#include <stan/math/prim/mat/prob/lkj_corr_cholesky_rng.hpp>
//#include <boost/random/mersenne_twister.hpp>

using Eigen::Dynamic;
using Eigen::Matrix;
using std::vector;

TEST(ProbDistributionsMultiNormalCholesky, MultiNormalVar) {
  using stan::math::var;
  Matrix<var, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<var, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<var, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  Matrix<var, Dynamic, Dynamic> L = Sigma.llt().matrixL();
  EXPECT_FLOAT_EQ(-11.73908,
                  stan::math::multi_normal_cholesky_log(y, mu, L).val());
}

TEST(AgradRev, check_varis_on_stack) {
  using stan::math::to_var;
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  Matrix<double, Dynamic, Dynamic> L = Sigma.llt().matrixL();
  test::check_varis_on_stack(stan::math::multi_normal_cholesky_log<true>(
      to_var(y), to_var(mu), to_var(L)));
  test::check_varis_on_stack(
      stan::math::multi_normal_cholesky_log<true>(to_var(y), to_var(mu), L));
  test::check_varis_on_stack(
      stan::math::multi_normal_cholesky_log<true>(to_var(y), mu, to_var(L)));
  test::check_varis_on_stack(
      stan::math::multi_normal_cholesky_log<true>(to_var(y), mu, L));
  test::check_varis_on_stack(
      stan::math::multi_normal_cholesky_log<true>(y, to_var(mu), to_var(L)));
  test::check_varis_on_stack(
      stan::math::multi_normal_cholesky_log<true>(y, to_var(mu), L));
  test::check_varis_on_stack(
      stan::math::multi_normal_cholesky_log<true>(y, mu, to_var(L)));

  test::check_varis_on_stack(stan::math::multi_normal_cholesky_log<false>(
      to_var(y), to_var(mu), to_var(L)));
  test::check_varis_on_stack(
      stan::math::multi_normal_cholesky_log<false>(to_var(y), to_var(mu), L));
  test::check_varis_on_stack(
      stan::math::multi_normal_cholesky_log<false>(to_var(y), mu, to_var(L)));
  test::check_varis_on_stack(
      stan::math::multi_normal_cholesky_log<false>(to_var(y), mu, L));
  test::check_varis_on_stack(
      stan::math::multi_normal_cholesky_log<false>(y, to_var(mu), to_var(L)));
  test::check_varis_on_stack(
      stan::math::multi_normal_cholesky_log<false>(y, to_var(mu), L));
  test::check_varis_on_stack(
      stan::math::multi_normal_cholesky_log<false>(y, mu, to_var(L)));
}

//  Here, we compare the speed of the new regression to that of one built from
//  existing primitives.
/*
TEST(ProbDistributionsMultiNormalCholesky, mvn_speed) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;

  typedef std::chrono::high_resolution_clock::time_point TimeVar;
#define duration(a) \
  std::chrono::duration_cast<std::chrono::microseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()

  // approx 2.3x faster for the new code
  const int R = 1500;
  const int C = 500;

  boost::random::mt19937 rng;

  int T1 = 0;
  int T2 = 0;
  for (size_t testnumber = 0; testnumber < 30; testnumber++) {
    std::vector<Matrix<double, Dynamic, 1>> Y_dbl;

    for (size_t i = 0; i < R; i++) {
      Matrix<double, Dynamic, 1> y_dbl = Eigen::VectorXd::Random(C);
      Y_dbl.push_back(y_dbl);
    }

    Matrix<double, Dynamic, 1> mu_dbl
        = Matrix<double, Dynamic, Dynamic>::Random(C, 1);

    Matrix<double, Dynamic, 1> sigma_dbl
        = Eigen::MatrixXd::Constant(C, 1, 0.01)
          + Matrix<double, Dynamic, 1>::Random(C, 1).array().abs().matrix();
    Matrix<double, Dynamic, Dynamic> L_Sigma_dbl
        = sigma_dbl.asDiagonal()
          * stan::math::lkj_corr_cholesky_rng(C, 1.0, rng);

    Matrix<var, Dynamic, 1> mu_v1 = mu_dbl;
    Matrix<var, Dynamic, Dynamic> L_Sigma_v1 = L_Sigma_dbl;
    std::vector<Matrix<var, Dynamic, 1>> Y_v1;
    for (size_t i = 0; i < R; i++) {
      Y_v1.push_back(Y_dbl[i]);
    }

    TimeVar t1 = timeNow();

    var lp1
        = stan::math::multi_normal_cholesky_old_lpdf(Y_v1, mu_v1, L_Sigma_v1);
    lp1.grad();

    TimeVar t2 = timeNow();
    stan::math::recover_memory();

    Matrix<var, Dynamic, 1> mu_v2 = mu_dbl;
    Matrix<var, Dynamic, Dynamic> L_Sigma_v2 = L_Sigma_dbl;
    std::vector<Matrix<var, Dynamic, 1>> Y_v2;
    for (size_t i = 0; i < R; i++) {
      Y_v2.push_back(Y_dbl[i]);
    }

    TimeVar t3 = timeNow();

    var lp2 = stan::math::multi_normal_cholesky_lpdf(Y_v2, mu_v2, L_Sigma_v2);
    lp2.grad();

    TimeVar t4 = timeNow();

    stan::math::recover_memory();

    T1 += duration(t2 - t1);
    T2 += duration(t4 - t3);
  }
  std::cout << "Existing Primitives:" << std::endl
            << T1 << std::endl
            << "New Primitives:" << std::endl
            << T2 << std::endl;
}
*/
