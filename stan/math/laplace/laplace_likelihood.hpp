#ifndef STAN_MATH_LAPLACE_LAPLACE_LIKELIHOOD_HPP
#define STAN_MATH_LAPLACE_LAPLACE_LIKELIHOOD_HPP

namespace stan {
namespace math {

// TO DO: create a parent structure, with each likelihood
// function acting as a child structure.
struct diff_poisson_log {
  Eigen::VectorXd n_samples_;
  Eigen::VectorXd sums_;

  diff_poisson_log();

  diff_poisson_log(const Eigen::VectorXd& n_samples,
                   const Eigen::VectorXd& sums)
    : n_samples_(n_samples), sums_(sums) { }

  template <typename T>
  void diff (const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta,
             Eigen::Matrix<T, Eigen::Dynamic, 1>& gradient,
             Eigen::Matrix<T, Eigen::Dynamic, 1>& hessian) const {
    hessian = - elt_multiply(n_samples_, exp(theta));
    gradient = sums_ + hessian;
  }
};

struct diff_logistic_log {
  Eigen::VectorXd n_samples_;
  Eigen::VectorXd sums_;

  diff_logistic_log();

  diff_logistic_log(const Eigen::VectorXd& n_samples,
                    const Eigen::VectorXd& sums)
    : n_samples_(n_samples), sums_(sums) { }

  template <typename T>
  T log_likelihood (const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta)
    const {
      Eigen::VectorXd one = rep_vector(1, theta.size());
      return sum(- sums_.cwiseProduct(log(one + exp(-theta)))
               - (n_samples_ - sums_).cwiseProduct(log(one + exp(theta))));

      // return dot_product(n_samples_, theta) - sum(log(exp(theta) + one));
    }

  template <typename T>
  void diff (const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta,
             Eigen::Matrix<T, Eigen::Dynamic, 1>& gradient,
             Eigen::Matrix<T, Eigen::Dynamic, 1>& hessian) const {
    gradient = sums_ - n_samples_.cwiseProduct(inv_logit(theta));

    Eigen::VectorXd exp_theta = exp(theta);
    Eigen::VectorXd one = rep_vector(1, theta.size());
    hessian = - n_samples_.cwiseProduct(elt_divide(exp_theta,
                                                   square(one + exp_theta)));

    // Eigen::VectorXd one = rep_vector(1, theta.size());
    // Eigen::VectorXd exp_theta = exp(theta);
    // Eigen::VectorXd one_plus_exp_theta = one + exp_theta;
    // Eigen::VectorXd one_plus_exp_neg_theta = one + exp(-theta);
    // 
    // gradient = sums_.cwiseProduct(elt_divide(one, one_plus_exp_theta)
    //   - (n_samples_ - sums_)
    //   .cwiseProduct(elt_divide(one, one_plus_exp_neg_theta)));
    // 
    // hessian = - n_samples_.cwiseProduct(elt_divide(exp_theta,
    //   square(one_plus_exp_theta)));

    // Eigen::VectorXd exp_theta = exp(theta);
    // Eigen::VectorXd one = rep_vector(1, theta.size());
    // hessian = - n_samples_.cwiseProduct(elt_divide(exp_theta,
    //                                                square(one + exp_theta)));
  }
  
  template <typename T>
  Eigen::Matrix<T, Eigen::Dynamic, 1>
  third_diff(const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta) const {
    Eigen::VectorXd exp_theta = exp(theta);
    Eigen::VectorXd one = rep_vector(1, theta.size());
    Eigen::VectorXd nominator = exp_theta.cwiseProduct(exp_theta - one);
    Eigen::VectorXd denominator = square(one + exp_theta)
                                   .cwiseProduct(one + exp_theta);

    return n_samples_.cwiseProduct(elt_divide(nominator, denominator));
  }
};

}  // namespace math
}  // namespace stan

#endif
