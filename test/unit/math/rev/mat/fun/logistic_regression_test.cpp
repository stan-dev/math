#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <vector>

struct lr {
  const std::vector<int>& y_;
  const Eigen::MatrixXd& x_;

  lr(const std::vector<int>& y, const Eigen::MatrixXd& x) : y_(y), x_(x) { }

  template <typename T>
  T operator()(const Eigen::Matrix<T, -1, 1>& beta) const {
    return stan::math::bernoulli_logit_regress_lpdf(y_, x_, beta);
  }
};

TEST(RevMat, LogisticRegressionOne) {
  using Eigen::VectorXd;
  using Eigen::MatrixXd;
  using Eigen::Matrix;
  using stan::math::var;
  using std::vector;

  int N = 100;
  int K = 10;

  vector<int> y(N, 1);
  VectorXd beta = VectorXd::Random(K);
  MatrixXd x = MatrixXd::Random(N, K);

  lr f(y, x);
  double val_ad = 0;
  VectorXd grad_ad;
  stan::math::gradient(f, beta, val_ad, grad_ad);


  double val_fd = 0;
  VectorXd grad_fd;
  stan::math::finite_diff_gradient(f, beta, val_fd, grad_fd, 1e-6);

  EXPECT_FLOAT_EQ(val_fd, val_ad);
  EXPECT_EQ(K, grad_fd.size());
  EXPECT_EQ(K, grad_ad.size());
  for (int k = 0; k < K; ++k)
    EXPECT_FLOAT_EQ(grad_fd(k), grad_ad(k));

}

