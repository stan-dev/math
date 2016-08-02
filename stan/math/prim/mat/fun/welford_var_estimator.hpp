#ifndef STANm_ATH_PRIMm_AT_FUN_WELFORD_VAR_ESTIMATOR_HPP
#define STANm_ATH_PRIMm_AT_FUN_WELFORD_VAR_ESTIMATOR_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <vector>

namespace stan {
  namespace math {

    class welford_var_estimator {
    public:
      explicit welford_var_estimator(int n)
        : m_(Eigen::VectorXd::Zero(n)),
          m2_(Eigen::VectorXd::Zero(n)) {
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
        m2_ += delta.cwiseProduct(q - m_);
      }

      int num_samples() { return _num_samples; }

      void sample_mean(Eigen::VectorXd& mean) { mean = m_; }

      void sample_variance(Eigen::VectorXd& var) {
        if (_num_samples > 1)
          var = m2_ / (_num_samples - 1.0);
      }

    protected:
      double _num_samples;

      Eigen::VectorXd m_;
      Eigen::VectorXd m2_;
    };

  }
}
#endif
