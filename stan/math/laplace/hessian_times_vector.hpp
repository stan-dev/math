#ifndef STAN_MATH_LAPLACE_HESSIAN_TIMES_VECTOR_HPP
#define STAN_MATH_LAPLACE_HESSIAN_TIMES_VECTOR_HPP

// TODO: refine include.
#include <stan/math/mix.hpp>

namespace stan {
namespace math {

  /**
   * Overload Hessian_times_vector function, under stan/math/mix/functor
   * to handle functions which take in arguments eta, delta, delta_int,
   * and pstream.
   */
  template <typename F>
  inline void hessian_times_vector(const F& f,
                            const Eigen::VectorXd& x,
                            const Eigen::VectorXd& eta,
                            const Eigen::VectorXd& delta,
                            const std::vector<int>& delta_int,
                            const Eigen::VectorXd& v,
                            double& fx,
                            Eigen::VectorXd& Hv,
                            std::ostream* pstream = 0) {
    using Eigen::Matrix;
    using Eigen::Dynamic;

    nested_rev_autodiff nested;

    int x_size = x.size();
    Matrix<var, Dynamic, 1> x_var = x;
    Matrix<fvar<var>, Dynamic, 1> x_fvar(x_size);
    for (int i = 0; i < x_size; i++) {
      x_fvar(i) = fvar<var>(x_var(i), v(i));
    }
    fvar<var> fx_fvar = f(x_fvar, eta, delta, delta_int, pstream);
    grad(fx_fvar.d_.vi_);
    Hv.resize(x_size);
    for (int i = 0; i < x_size; i++) Hv(i) = x_var(i).adj();
}

}  // namespace math
}  // namespace stan

#endif
