#ifndef STAN_MATH_REV_MAT_FUN_MATRIX_EXP_ACTION_HANDLER_HPP
#define STAN_MATH_REV_MAT_FUN_MATRIX_EXP_ACTION_HANDLER_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <vector>

namespace stan {
namespace math {

  double l1norm(const Eigen::MatrixXd& m) {
    return m.colwise().lpNorm<1>().maxCoeff();
  }

  // double matrix_infty_norm(const Eigen::MatrixXd& m) {
  //   return m.rowwise().lpNorm<1>().maxCoeff();
  // }

  class matrix_exp_action_handler {
    static constexpr int p_max = 8;
    static constexpr int m_max = 55;
    static constexpr double tol = 1.1e-16;
    static const std::vector<double> theta_m_single_precision;
    static const std::vector<double> theta_m_double_precision;

  public:
    matrix_exp_action_handler() {}

    inline
    Eigen::MatrixXd action(const Eigen::MatrixXd& mat,
                           const Eigen::MatrixXd& b,
                           const double& t = 1.0) {
      Eigen::MatrixXd A = mat;
      double mu = A.trace()/A.rows();
      for (int i = 0; i < A.rows(); ++i) {
        A(i, i) -= mu;
      }

      int m{0}, s{0};
      set_approximation_parameter(A, t, m, s);

      Eigen::MatrixXd res(A.rows(), b.cols());

      for (int col = 0; col < b.cols(); ++col) {
        bool conv = false;
        Eigen::VectorXd B = b.col(col);
        Eigen::VectorXd F = B;
        const auto eta = std::exp(t * mu / s);
        for (int i = 1; i < s+1; ++i) {
          auto c1 = B.template lpNorm<Eigen::Infinity>();
          if (m > 0) {
            for (int j=1; j < m+1; ++j) {
              B = t * A * B / (s * j);
              auto c2 = B.template lpNorm<Eigen::Infinity>();
              F += B;
              if (c1 + c2 < tol * F.template lpNorm<Eigen::Infinity>()) {
                conv = true;
                break;
              }
              c1 = c2;
            }
          }
          F *= eta;
          B = F;
          if (conv) break;
        }
        res.col(col) = F;
      }     // loop b columns
      return res;
    }

    inline
    void
    set_approximation_parameter(const Eigen::MatrixXd& mat,
                                const double& t, int& m, int& s) {
      // to simply the process, always use m_max
      if (l1norm(mat) < tol || t < tol) {
        m = 0;
        s = 1;
      } else {
        m = m_max;

        Eigen::MatrixXd a = mat * t;
        const double theta_m = theta_m_double_precision.back();
        for (auto i = 0; i < std::ceil(std::log2(p_max)); ++i) {
          a *= a;
        }
        double ap = std::pow(l1norm(a), 1.0/p_max);
        int c = std::ceil(ap/theta_m);
        s = (c < 1 ? 1 : c);
      }
    }
  };

  const std::vector<double>
  matrix_exp_action_handler::theta_m_single_precision{
    1.3e-1, 1.0e0, 2.2e0, 3.6e0, 4.9e0, 6.3e0, 7.7e0, 9.1e0,
      1.1e1, 1.2e1, 1.3e1};

  const std::vector<double>
  matrix_exp_action_handler::theta_m_double_precision{
    2.4e-3, 1.4e-1, 6.4e-1, 1.4e0, 2.4e0, 3.5e0, 4.7e0,
      6.0e0, 7.2e0, 8.5e0, 9.9e0};

}
}

#endif
