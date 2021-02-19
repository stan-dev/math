#include <stan/math.hpp>
#include <stan/math/laplace/laplace_likelihood.hpp>
#include <stan/math/laplace/third_diff_directional.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/mix.hpp>
#include <stan/math/fwd.hpp>
// #include <stan/math/mix/functor/hessian_times_vector.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>

// This is what a function define in Stan would return.
struct neg_bin_log_likelihood {
  template <typename T_theta, typename T_eta>
  stan::return_type_t<T_theta, T_eta>
  operator()(const Eigen::Matrix<T_theta, -1, 1>& theta,
             const Eigen::Matrix<T_eta, -1, 1>& eta,
             const Eigen::VectorXd& delta,
             const std::vector<int>& delta_int,
             std::ostream* pstream) const {
    stan::math::diff_neg_binomial_2_log
      diff_functor(delta, delta_int, theta.size());

    return diff_functor.log_likelihood(theta, eta);
  }
};

template <typename F>
struct f_theta {
  Eigen::VectorXd eta_;
  Eigen::VectorXd delta_;
  std::vector<int> delta_int_;
  std::ostream* pstream_;
  F f_functor_;

  f_theta(const Eigen::VectorXd& eta,
          const Eigen::VectorXd& delta,
          const std::vector<int>& delta_int,
          std::ostream* pstream,
          F f_functor) :
    eta_(eta), delta_(delta), delta_int_(delta_int), f_functor_(f_functor) { }

  template <typename T>
  T operator()(const Eigen::Matrix<T, -1, 1>& theta) const {
    return f_functor_(theta, eta_, delta_, delta_int_, pstream_);
  }
};

TEST(laplace_diff, gradient) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::value_of;
  using stan::math::hessian_times_vector;
  using stan::math::nested_rev_autodiff;
  using stan::math::third_diff_directional;

  Eigen::Matrix<var, -1, 1> theta(2);
  theta << 1, 1;
  Eigen::VectorXd theta_dbl = value_of(theta);
  Eigen::VectorXd eta(1);
  eta << 1.2;

  Eigen::VectorXd y(2);
  y << 1, 1;
  std::vector<int> y_index(2);
  y_index[0] = 0;
  y_index[1] = 1;

  neg_bin_log_likelihood likelihood;
  f_theta<neg_bin_log_likelihood> f(eta, y, y_index, 0, likelihood);

  // var log_density = likelihood(theta, eta, y, y_index, 0);
  var log_density = f(theta);
  std::cout << log_density.val() << std::endl;

  // autodiff for first derivative
  // if (FALSE) {
  //   VEC g;
  //   AVEC parm_vec = createAVEC(theta(0), theta(1));
  //   log_density.grad(parm_vec, g);
  //   std::cout << "Gradien: ";
  //   for (int i = 0; i < theta.size(); i++) std::cout << g[i] << " ";
  //   std::cout << std::endl;
  // }

  // specify initial tangent
  Eigen::VectorXd tangent(2);
  tangent << 1, 1;

  double fx;
  Eigen::VectorXd Hv;

  hessian_times_vector(f, theta_dbl, tangent, fx, Hv);

  std::cout << "value: " << fx << std::endl;
  std::cout << "hessian-vector: " << Hv.transpose() << std::endl;

  // Compute third-order derivative
  if (TRUE) {
    using Eigen::Matrix;
    nested_rev_autodiff nested;

    // CHECK -- why would need a for loop for assignment (see
    // hessian_times_vector.hpp).
    Matrix<var, -1, 1> theta_var = theta_dbl;

    Matrix<fvar<var>, -1, 1> theta_fvar(theta_dbl.size());
    for (int i = 0; i < theta_dbl.size(); ++i) {
      theta_fvar(i) = fvar<var>(theta_var(i), tangent(i));
    }
    fvar<var> fx_fvar = f(theta_fvar);

    std::cout << "fx: " << value_of(fx_fvar.val_) << std::endl;
    std::cout << "grad_fx_dot_v: " << value_of(fx_fvar.d_) << std::endl;

    // grad(fx_fvar.d_.vi_);
    // for (int i = 0; i < theta_dbl.size(); ++i)
    //   std::cout << theta_var(i).adj() << " ";
    // std::cout << std::endl;

    Matrix<fvar<fvar<var>>, -1, 1> theta_ffvar(theta_dbl.size());
    for (int i = 0; i < theta_dbl.size(); ++i) {
      theta_ffvar(i) = fvar<fvar<var>>(theta_fvar(i), tangent(i));
    }
    fvar<fvar<var>> fx_ffvar = f(theta_ffvar);
    var grad2_fx_dot_vv = fx_ffvar.d_.d_;

    std::cout << "fx: " << value_of(fx_ffvar.val_.val_) << std::endl;
    std::cout << "grad_grad_fx_dot_v: "
      << value_of(fx_ffvar.d_.d_) << std::endl;

    // var fx_var;
    // var grad3_f_v;
    grad(grad2_fx_dot_vv.vi_);

    Eigen::VectorXd grad3_f_v(theta_dbl.size());
    for (int i = 0; i < theta_dbl.size(); ++i)
      grad3_f_v(i) = theta_var(i).adj();

    std::cout << "grad3 f: " << grad3_f_v.transpose() << std::endl;

    // var fx_var;
    // var grad_fx_var_dot_v;
    // gradient_dot_vector(f, theta_var, tang, fx_var, grad_fx_var_dot_v);
    // fx = fx_var.val();
    // grad(grad_fx_var_dot_v.vi_);
    // Hv.resize(theta.size());
    // for (int i = 0; i < theta.size(); ++i) Hv(i) = theta_var(i).adj();
  }

  // Test function
  Eigen::VectorXd third_diff;
  third_diff_directional(f, theta_dbl, fx, third_diff, tangent, tangent);

  std::cout << "f: " << fx << std::endl;
  std::cout << "third diff: " << third_diff.transpose() << std::endl;

}
