#ifndef STANm_ATH_PRIMm_AT_FUN_WELFORD_COVAR_ESTIMATOR_HPP
#define STANm_ATH_PRIMm_AT_FUN_WELFORD_COVAR_ESTIMATOR_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <vector>

namespace stan {
  namespace math {

    class welford_covar_estimator {
    public:
      explicit welford_covar_estimator(int n)
        : m_(Eigen::VectorXd::Zero(n)),
          m2_(Eigen::MatrixXd::Zero(n, n)) {
        restart();
      }

      void restart() {
        _num_samples = 0;
        m_.setZero();
        m2_.setZero();
      }

      void add_sample(const Eigen::VectorXd& q) {
        ++_num_samples;

        Eigen::VectorXd delta(q - m_);
        m_  += delta / _num_samples;
        m2_ += (q - m_) * delta.transpose();
      }

      int num_samples() { return _num_samples; }

      void sample_mean(Eigen::VectorXd& mean) { mean = m_; }

      void sample_covariance(Eigen::MatrixXd& covar) {
        if (_num_samples > 1)
          covar = m2_ / (_num_samples - 1.0);
      }

    protected:
      double _num_samples;

      Eigen::VectorXd m_;
      Eigen::MatrixXd m2_;
    };

  }
}
#endif
