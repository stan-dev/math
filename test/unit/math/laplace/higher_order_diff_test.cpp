#include <stan/math.hpp>
#include <stan/math/laplace/laplace_likelihood.hpp>
#include <stan/math/laplace/laplace_likelihood_general.hpp>
#include <stan/math/laplace/third_diff_directional.hpp>
#include <stan/math/laplace/hessian_times_vector.hpp>
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

  f_theta() { }  // default constructor required for default class.

  f_theta(const F& f_functor,
          const Eigen::VectorXd& eta,
          const Eigen::VectorXd& delta,
          const std::vector<int>& delta_int,
          std::ostream* pstream) :
    eta_(eta), delta_(delta), delta_int_(delta_int), f_functor_(f_functor) { }

  template <typename T>
  T operator()(const Eigen::Matrix<T, -1, 1>& theta) const {
    return f_functor_(theta, eta_, delta_, delta_int_, pstream_);
  }
};

class neg_bin_log_diff_test : public::testing::Test {
protected:
  void SetUp() override {
    theta.resize(2);
    theta << 1, 1;
    theta_dbl = value_of(theta);
    eta.resize(1);
    eta << 1.2;

    y.resize(2);
    y << 0, 1;
    y_index.resize(2);
    y_index[0] = 0;
    y_index[1] = 1;

    f_theta<neg_bin_log_likelihood> f_(likelihood, eta, y, y_index, 0);
    f = f_;
  }

  Eigen::Matrix<stan::math::var, -1, 1> theta;
  Eigen::VectorXd theta_dbl;
  Eigen::VectorXd eta;
  Eigen::VectorXd y;
  std::vector<int> y_index;
  neg_bin_log_likelihood likelihood;
  f_theta<neg_bin_log_likelihood> f;
};


TEST_F(neg_bin_log_diff_test, manual_calls) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::value_of;
  using stan::math::hessian_times_vector;
  using stan::math::nested_rev_autodiff;
  using stan::math::third_diff_directional;

  // var log_density = likelihood(theta, eta, y, y_index, 0);
  var log_density = f(theta);
  std::cout << log_density.val() << std::endl;

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
  }

  // Test function for directional Hessian and directional third diff.
  Eigen::VectorXd hessian_v;
  hessian_times_vector(likelihood, theta_dbl, eta, y, y_index,
                       tangent, fx, hessian_v, 0);

  std::cout << "hessian_v: " << hessian_v.transpose() << std::endl;

  Eigen::VectorXd third_diff;
  third_diff_directional(likelihood, theta_dbl, eta, y, y_index,
                         fx, third_diff, tangent, tangent, 0);

  std::cout << "f: " << fx << std::endl;
  std::cout << "third diff: " << third_diff.transpose() << std::endl;
}

TEST_F(neg_bin_log_diff_test, diff_likelihood) {
  using stan::math::diff_likelihood;
  using Eigen::VectorXd;

  diff_likelihood<neg_bin_log_likelihood> lk(likelihood, y, y_index, 0);
  double lpmf = lk.log_likelihood(theta_dbl, eta);
  VectorXd gradient, hessian;
  lk.diff(theta_dbl, eta, gradient, hessian);
  VectorXd third_diff = lk.third_diff(theta_dbl, eta);

  std::cout << "lpmf: " << lpmf << std::endl
    << "gradient: " << gradient.transpose() << std::endl
    << "hessian: " << hessian.transpose() << std::endl
    << "third diff: " << third_diff.transpose() << std::endl;

}
