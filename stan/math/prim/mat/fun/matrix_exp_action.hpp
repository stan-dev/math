#ifndef STAN_MATH_PRIM_MAT_FUN_MATRIX_EXP_ACTION_HPP
#define STAN_MATH_PRIM_MAT_FUN_MATRIX_EXP_ACTION_HPP

namespace stan {
namespace math {

  template <typename T>
  double matrix_1_norm(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m) {
    auto m_dbl {to_var(m)};
    return matrix_1_norm(m_dbl);
  }

  template <>
  double matrix_1_norm(const Eigen::MatrixXd& m) {
    return m.colwise().lpNorm<1>().maxCoeff();
  }

  template <typename T>
  double matrix_infty_norm(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m) {
    auto m_dbl {to_var(m)};
    return matrix_infty_norm(m_dbl);
  }

  template <>
  double matrix_infty_norm(const Eigen::MatrixXd& m) {
    return m.rowwise().lpNorm<1>().maxCoeff();
  }

  class matrix_exp_action_handler {

    static constexpr int p_max = 8;
    static constexpr int m_max = 55;
    static constexpr double tol = 1.1e-16;
    static const std::vector<double> theta_m_single_precision;
    static const std::vector<double> theta_m_double_precision;

    void set_approximation_parameter(const Eigen::MatrixXd& mat, int& m, int& s) {
      // to simply the process, always use m_max
      if (matrix_1_norm(mat) < tol) {
        m = 0;
        s = 1;
      } else {
        m = m_max;

        Eigen::MatrixXd a = mat;
        const double theta_m = theta_m_double_precision.back();
        int c = std::ceil(matrix_1_norm(a)/theta_m);
        for (size_t i = 1; i < p_max; ++i) {
          a *= mat;
          double ap = std::pow(matrix_1_norm(a), 1.0/i);
          int ca = std::ceil(ap/theta_m);
          if (ca < c) {
            c = ca;
          }
        }
        s = (c < 1 ? 1 : c);        
      }
    }

    Eigen::VectorXd action(const Eigen::MatrixXd& mat, const
                     Eigen::VectorXd& b, const double& t) {
      Eigen::MatrixXd A = mat;
      Eigen::VectorXd B = b;
      double mu = A.trace()/A.rows();
      for (int i = 0; i < A.rows(); ++i) {
        A(i, i) -= mu;
      }

      int m{0}, s{0};
      set_approximation_parameter(A, m, s);
      auto F = B;
      const auto eta = std::exp(t * mu / s);
      for (int i = 1; i < s+1; ++i) {
        auto c1 = B.lpNorm<Eigen::Infinity>();
        if (m > 0) {
          for (int j=1; j < m+1; ++j) {
            B = t * A * B / (s * j);
            auto c2 = B.lpNorm<Eigen::Infinity>();
            F += B;
            if (c1 + c2 < tol * F.lpNorm<Eigen::Infinity>()) {
              break;
            }
            c1 = c2;
          }
        }
        F *= eta;
        B = F;          
      }
      return F;
    }

  public:
    matrix_exp_action_handler() {}
    Eigen::VectorXd perform_action(const Eigen::MatrixXd& mat, const
                                   Eigen::VectorXd& b, const double& t) {
      return action(mat, b, t);
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

}  // namespace math
}  // namespace stan
#endif
