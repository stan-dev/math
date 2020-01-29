#ifndef STAN_MATH_LAPLACE_LAPLACE_PSEUDO_TARGET_HPP
#define STAN_MATH_LAPLACE_LAPLACE_PSEUDO_TARGET_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/quad_form.hpp>
#include <iostream>

namespace stan {
namespace math {
  /**
   * Function to compute the pseudo target, $\tilde Z$,
   * with a custom derivative method.
   */
   inline double laplace_pseudo_target (
     const Eigen::MatrixXd& K,
     const Eigen::VectorXd& a,
     const Eigen::MatrixXd& R,
     const Eigen::VectorXd& l_grad,
     const Eigen::VectorXd& s2) {
       double s1 = 0.5 * quad_form(K, a) - 0.5 * sum((R * K).diagonal());
       Eigen::VectorXd b = K * l_grad;
       Eigen::VectorXd s3 = b - K * (R * b);
       return s1 + s2.dot(s3);
     }

  /**
   * Vari class for the function.
   */
  struct laplace_pseudo_target_vari : public vari {
    /* number of elements in covariance matrix. */
    int K_size_;
    /* covariance matrix. */
    vari** K_;
    /* pseudo target. */
    vari** pseudo_target_;
    /* An object to store the sensitivities of K. */
    Eigen::MatrixXd K_adj_;

    template <typename T>
    laplace_pseudo_target_vari (
      const Eigen::VectorXd& a,
      const Eigen::MatrixXd& R,
      const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& K,
      const Eigen::VectorXd& s2,
      const Eigen::VectorXd& l,
      double pseudo_target)
      : vari(pseudo_target),
        K_size_(K.size()),
        K_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(
          K.size())),
        pseudo_target_(
          ChainableStack::instance_->memalloc_.alloc_array<vari*>(1)) {
        int dim_theta = K.rows();
        for (int j = 0; j < dim_theta; j++)
          for (int i = 0; i < dim_theta; i++)
            K_[j * dim_theta + i] = K(i, j).vi_;

        pseudo_target_[0] = this;
        pseudo_target_[0] = new vari(pseudo_target, false);

        K_adj_ = 0.5 * a * a.transpose() - 0.5 * R
           + s2 * l.transpose()
           - (R * (value_of(K) * s2)) * l.transpose();
      }

      void chain() {
        int dim_theta = K_adj_.rows();
        for (int j = 0; j < dim_theta; j++)
           for (int i = 0; i < dim_theta; i++)
             K_[j * dim_theta + i]->adj_ +=
               pseudo_target_[0]->adj_ * K_adj_(i, j);
     }
  };

  /**
   * Overload function for case where K is passed as a matrix of var.
   */
  template <typename T>
    inline T laplace_pseudo_target (
      const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& K,
      const Eigen::VectorXd& a,
      const Eigen::MatrixXd& R,
      const Eigen::VectorXd& l_grad,
      const Eigen::VectorXd& s2) {
      double pseudo_target_dbl
        = laplace_pseudo_target(value_of(K), a, R, l_grad, s2);

      // construct vari
      laplace_pseudo_target_vari* vi0
        = new laplace_pseudo_target_vari(a, R, K, s2, l_grad,
                                         pseudo_target_dbl);

      var pseudo_target = var(vi0->pseudo_target_[0]);
      return pseudo_target;
  }

}  // namespace math
}  // namespace stan

#endif
