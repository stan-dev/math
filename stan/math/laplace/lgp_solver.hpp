#ifndef STAN_MATH_LAPLACE_LGP_SOLVER_HPP
#define STAN_MATH_LAPLACE_LGP_SOLVER_HPP

#include <stan/math/laplace/lgp_dense_system.hpp>


namespace stan {
namespace math {

/**
 * A functor that returns f, the system function, for
 * a latent Gaussian Poisson process. Note that the
 * number of samples per group, sum per group,
 * and Gaussian precision matrix (Q) are pre-computed.
 */
struct lgp_f {
  inline Eigen::VectorXd
  operator() (const Eigen::VectorXd& theta,
              const Eigen::VectorXd& phi,
              const std::vector<double>& dat,
              const std::vector<int>& dat_int,
              std::ostream* pstream__) const {
    size_t N = theta.size();
    Eigen::Map<const Eigen::VectorXd> n_samples(&dat[0], N);
    Eigen::Map<const Eigen::VectorXd> sums(&dat[N], N);
    Eigen::Map<const Eigen::MatrixXd> Q(&dat[2 * N], N, N);

    return sums - n_samples.cwiseProduct(exp(theta)) - Q * theta;
  }
};

/**
 * A functor that returns the Jacobian of f, the system function,
 * with respect to the unknown theta.
 * CHECK -- do I need that strict a signature, or can I relax it?
 * FIX ME -- why do we need templates for x to make elt_multiply work?
 */
struct lgp_J_f {
  template <typename F>
  inline int
  operator() (const F& f,
              const Eigen::VectorXd& theta,
              const Eigen::VectorXd& phi,
              const std::vector<double>& dat,
              const std::vector<int>& dat_int,
              std::ostream* msgs,
              const double x_sun[], SUNMatrix J) const {
    size_t N = theta.size();
    // const std::vector<double> x_vec(x_sun, x_sun + N);

    Eigen::Map<const Eigen::VectorXd> n_samples(&dat[0], N);
    Eigen::Map<const Eigen::MatrixXd> Q(&dat[2 * N], N, N);
    Eigen::Map<const Eigen::VectorXd> x(&x_sun[0], N);

    std::vector<double> jacobian_x = std::vector<double>(N * N);
    Eigen::Map<Eigen::MatrixXd>(&jacobian_x[0], N, N)
      = (- n_samples).cwiseProduct(exp(x)).asDiagonal();
    Eigen::Map<Eigen::MatrixXd>(&jacobian_x[0], N, N) -= Q;

    std::move(jacobian_x.begin(), jacobian_x.end(), SM_DATA_D(J));

    return 0;
  }
};


}  // namespace math
}  // namespace stan

#endif
