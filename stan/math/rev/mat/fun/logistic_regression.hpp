#ifndef STAN_MATH_REV_MAT_FUN_LOGISTIC_REGRESSION_HPP
#define STAN_MATH_REV_MAT_FUN_LOGISTIC_REGRESSION_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/scal/fun/inv_logit.hpp>
#include <stan/math/prim/scal/fun/log_inv_logit.hpp>
#include <stan/math/prim/scal/fun/log1m.hpp>
#include <stan/math/prim/scal/fun/log1m_inv_logit.hpp>
#include <vector>
#include <cmath>

namespace stan {
  namespace math {

    // usable even in a hierarchical logistic regression setting

    // BETTER IF WE CAN SORT OUT 0/1 outcomes in tdata
    // logistic_regression(matrix x1, matrix x0, vector beta)

    //  partials = ((y- y_hat).asDiagonal() * x).colwise().sum();

    //  partials = ((-y_hat).asDiagonal() * x1).colwise().sum()
    //             + ((1 - y_hat).asDiagonal() * x0).colwise().sum();

    //  partials = (-y_hat.asDiagonal() * x).colwise().sum()
    //             + x0.colwise().sum();




    // MORE STABLE WITH:
    // log_inv_logit and log1m_inv_logit in val calc

    // SIMPLER WITH VECTORIZED inv_logit

    // y :    N x 1
    // x:     N x K
    // beta:  K x 1    x * beta: N x 1
    template <bool propto>
    var bernoulli_logit_regress_lpdf(const std::vector<int>& y,
                            const Eigen::MatrixXd& x,
                            const Eigen::Matrix<var, -1, 1>& beta) {
      using Eigen::VectorXd;
      using Eigen::RowVectorXd;

      int N = x.rows();
      int K = x.cols();

      Eigen::VectorXd beta_d(K);
      for (int k = 0; k < K; ++k)
        beta_d(k) = beta(k).val();

      VectorXd logit_y_hat = x * beta_d;  // mat multiply to speedup
      VectorXd y_hat(N);
      for (int n = 0; n < N; ++n)
        y_hat(n) = inv_logit(logit_y_hat(n));

      double val = 0;
      for (int n = 0; n < N; ++n)
        val += y[n] ? std::log(y_hat[n]) : log1m(y_hat[n]);

      VectorXd y_d(N);
      for (int n = 0; n < N; ++n) y_d(n) = y[n];
      VectorXd partials = ((y_d - y_hat).asDiagonal() * x).colwise().sum();
      // RowVectorXd partials = RowVectorXd::Zero(K);
      // for (int n = 0; n < N; ++n)
      // partials += x.row(n) * (y[n] - y_hat[n]);

      vari** operands = ChainableStack::memalloc_.alloc_array<vari*>(K);
      for (int k = 0; k < K; ++k)
        operands[k] = beta(k).vi_;

      double* parts = ChainableStack::memalloc_.alloc_array<double>(K);
      for (int k = 0; k < K; ++k)
        parts[k] = partials[k];  // TODO(carpenter): use parts directly

      return var(new stored_gradient_vari(val, K, operands, parts));
    }

    template <bool propto>
    double bernoulli_logit_regress_lpdf(const std::vector<int>& y,
                                        const Eigen::MatrixXd& x,
                                        const Eigen::VectorXd& beta_d) {
      using Eigen::VectorXd;
      int N = x.rows();
      VectorXd logit_y_hat = x * beta_d;  // mat multiply to speedup
      VectorXd y_hat(N);
      for (int n = 0; n < N; ++n)
        y_hat(n) = inv_logit(logit_y_hat(n));
      double val = 0;
      for (int n = 0; n < N; ++n)
        val += y[n]
          ? log_inv_logit(logit_y_hat[n])
          : log1m_inv_logit(logit_y_hat[n]);
      return val;
    }

  }
}
#endif
