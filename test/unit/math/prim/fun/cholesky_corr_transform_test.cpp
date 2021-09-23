#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(ProbTransform, CholeskyCorrelation4) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, Dynamic> L(4, 4);
  L << 1, 0, 0, 0, -0.2, 0.9797959, 0, 0, 0.5, -0.3, 0.8124038, 0, 0.7, -0.2,
      0.6, 0.3316625;

  Matrix<double, Dynamic, 1> y = stan::math::cholesky_corr_free(L);

  Matrix<double, Dynamic, Dynamic> x
      = stan::math::cholesky_corr_constrain(y, 4);

  Matrix<double, Dynamic, 1> yrt = stan::math::cholesky_corr_free(x);

  EXPECT_EQ(y.size(), yrt.size());
  for (int i = 0; i < yrt.size(); ++i)
    EXPECT_FLOAT_EQ(y(i), yrt(i));

  for (int m = 0; m < 4; ++m)
    for (int n = 0; n < 4; ++n)
      EXPECT_FLOAT_EQ(L(m, n), x(m, n));
}

void test_cholesky_correlation_values(
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& L) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::cholesky_corr_constrain;
  using stan::math::cholesky_corr_free;
  using std::vector;
  double lp = 0;
  int K = L.rows();
  int K_choose_2 = (K * (K - 1)) / 2;

  // test number of free parameters
  Matrix<double, Dynamic, 1> y = stan::math::cholesky_corr_free(L);
  EXPECT_EQ(K_choose_2, y.size());

  std::vector<Matrix<double, Dynamic, 1>> y_vec
      = stan::math::cholesky_corr_free(
          std::vector<Matrix<double, Dynamic, Dynamic>>{L, L, L});
  for (int i = 0; i < y_vec.size(); ++i) {
    EXPECT_EQ(K_choose_2, y_vec[i].size());
  }

  // test transform roundtrip without Jacobian
  std::vector<Matrix<double, Dynamic, Dynamic>> x_vec
      = stan::math::cholesky_corr_constrain<false>(y_vec, K, lp);

  std::vector<Matrix<double, Dynamic, 1>> yrt
      = stan::math::cholesky_corr_free(x_vec);

  EXPECT_EQ(y.size(), yrt[0].size());
  for (int i = 0; i < yrt.size(); ++i) {
    for (int j = 0; j < yrt[i].size(); ++j) {
      EXPECT_FLOAT_EQ(y_vec[i](j), yrt[i](j));
    }
  }
  for (auto&& x_i : x_vec) {
    for (int m = 0; m < K; ++m) {
      for (int n = 0; n < K; ++n) {
        EXPECT_FLOAT_EQ(L(m, n), x_i(m, n));
      }
    }
  }

  // test transform roundtrip with Jacobian (Jacobian itself tested above)
  std::vector<Matrix<double, Dynamic, Dynamic>> x2_vec
      = stan::math::cholesky_corr_constrain<true>(y_vec, K, lp);

  std::vector<Matrix<double, Dynamic, 1>> yrt2
      = stan::math::cholesky_corr_free(x2_vec);

  for (auto&& yrt_i : yrt2) {
    EXPECT_EQ(y.size(), yrt_i.size());
    for (int i = 0; i < yrt_i.size(); ++i) {
      EXPECT_FLOAT_EQ(y(i), yrt_i(i));
    }
  }
  for (auto&& x2_i : x2_vec) {
    for (int m = 0; m < K; ++m) {
      for (int n = 0; n < K; ++n) {
        EXPECT_FLOAT_EQ(L(m, n), x2_i(m, n));
      }
    }
  }
}

TEST(ProbTransform, CholeskyCorrelationRoundTrips) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  Matrix<double, Dynamic, Dynamic> L1(1, 1);
  L1 << 1;
  test_cholesky_correlation_values(L1);

  Matrix<double, Dynamic, Dynamic> L2(2, 2);
  L2 << 1, 0, -0.5, 0.8660254;
  test_cholesky_correlation_values(L2);

  Matrix<double, Dynamic, Dynamic> L4(4, 4);
  L4 << 1, 0, 0, 0, -0.2, 0.9797959, 0, 0, 0.5, -0.3, 0.8124038, 0, 0.7, -0.2,
      0.6, 0.3316625;
  test_cholesky_correlation_values(L4);
}
