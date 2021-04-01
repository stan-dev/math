#include <stan/math.hpp>
#include <stan/math/laplace/laplace.hpp>
#include <stan/math/laplace/laplace_likelihood_general.hpp>
#include <stan/math/laplace/laplace_likelihood_poisson_log.hpp>
#include <stan/math/laplace/prob/laplace_base_rng.hpp>
#include <stan/math/laplace/prob/laplace_poisson_log_rng.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>

#include <test/unit/math/laplace/laplace_utility.hpp>
#include <test/unit/math/rev/fun/util.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>

struct poisson_log_likelihood {
  template<typename T_theta, typename T_eta>
  stan::return_type_t<T_theta, T_eta>
  operator()(const Eigen::Matrix<T_theta, -1, 1>& theta,
             const Eigen::Matrix<T_eta, -1, 1>& eta,
             const Eigen::VectorXd& delta,
             const std::vector<int>& n_samples,
             std::ostream* pstream) const {
    using stan::math::to_vector;
    using stan::math::log;
    int n = 911;
    Eigen::VectorXd y = delta.head(n);
    Eigen::VectorXd ye = delta.tail(n);
    // Eigen::VectorXd log_ye = ye.log();

    stan::math::diff_poisson_log
      diff_functor(to_vector(n_samples), y, log(ye));

    return diff_functor.log_likelihood(theta, eta);
  }
};

// TODO(charlesm93): update using new function signatures.
class laplace_disease_map_test : public::testing::Test {
protected:
  void SetUp() override {
    dim_theta = 911;
    n_observations = 911;
    data_directory = "test/unit/math/laplace/aki_disease_data/";
    x1.resize(dim_theta);
    x2.resize(dim_theta);
    y.resize(n_observations);
    ye.resize(n_observations);
    read_in_data(dim_theta, n_observations, data_directory, x1, x2, y, ye);

    if (FALSE) {
      // look at some of the data
      std::cout << "x_1: " << x1[0] << " " << x2[0] << std::endl
                << "x_2: " << x1[1] << " " << x2[1] << std::endl
                << "y_1: " << y[0] << " y_2: " << y[1] << std::endl
                << "ye_1: " << ye[0] << " ye_2: " << ye[1] << std::endl;
    }

    dim_x = 2;
    x.resize(dim_theta);
    for (int i = 0; i < dim_theta; i++) {
      Eigen::VectorXd coordinate(dim_x);
      coordinate << x1[i], x2[i];
      x[i] = coordinate;
    }

    // one observation per group
    n_samples.resize(dim_theta);
    for (int i = 0; i < dim_theta; i++) n_samples[i] = 1;

    theta_0 = Eigen::VectorXd::Zero(dim_theta);
    dim_phi = 2;
    phi.resize(dim_phi);
    phi << 0.3162278, 200;  // variance, length scale

    delta_lk.resize(2 * n_observations);
    for (int i = 0; i < n_observations; i++) delta_lk(i) = y[i];
    for (int i = 0; i < n_observations; i++)
      delta_lk(n_observations + i) = ye(i);
  }

  int dim_theta;
  int n_observations;
  std::string data_directory;
  std::vector<double> x1, x2;
  std::vector<int> y;
  Eigen::VectorXd ye;
  int dim_x;
  std::vector<Eigen::VectorXd> x;
  std::vector<int> n_samples;
  std::vector<double> delta;
  std::vector<int> delta_int;

  Eigen::VectorXd theta_0;
  int dim_phi;
  Eigen::Matrix<stan::math::var, -1, 1> phi;
  Eigen::Matrix<stan::math::var, -1, 1> eta_dummy;

  Eigen::VectorXd delta_lk;
  poisson_log_likelihood f;
};


TEST_F(laplace_disease_map_test, lk_analytical) {

  // Based on (Vanhatalo, Pietilainen and Vethari, 2010). See
  // https://research.cs.aalto.fi/pml/software/gpstuff/demo_spatial1.shtml
  using stan::math::var;
  using stan::math::laplace_marginal_poisson_log_lpmf;

  auto start = std::chrono::system_clock::now();

  var marginal_density
    = laplace_marginal_poisson_log_lpmf(y, n_samples, ye,
                                        sqr_exp_kernel_functor(),
                                        phi, x, delta, delta_int, theta_0);

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = end - start;

  VEC g;
  AVEC parm_vec = createAVEC(phi(0), phi(1));
  marginal_density.grad(parm_vec, g);

  std::cout << "LAPLACE MARGINAL AND VARI CLASS" << std::endl
            << "density: " << value_of(marginal_density) << std::endl
            << "autodiff grad: " << g[0] << " " << g[1] << std::endl
            << "total time: " << elapsed_time.count() << std::endl
            << std::endl;

  // Expected result
  // density: -2866.88
  // autodiff grad: 266.501 -0.425901
  // total time: 0.414122 (on new computer), 0.627501 (on old computer)

  // TODO(charlesm93): update signatures for rng functions.
  ////////////////////////////////////////////////////////////////////////
  // Let's now generate a sample theta from the estimated posterior

  using stan::math::diff_poisson_log;
  using stan::math::to_vector;
  using stan::math::laplace_base_rng;
  using stan::math::laplace_poisson_log_rng;

  diff_poisson_log diff_likelihood(to_vector(n_samples),
                                   to_vector(y),
                                   stan::math::log(ye));
  boost::random::mt19937 rng;
  start = std::chrono::system_clock::now();
  Eigen::VectorXd
    theta_pred = laplace_base_rng(diff_likelihood,
                                  sqr_exp_kernel_functor(),
                                  phi, eta_dummy, x, x, delta, delta_int,
                                  theta_0, rng);

  end = std::chrono::system_clock::now();
  elapsed_time = end - start;

  std::cout << "LAPLACE_APPROX_RNG" << std::endl
            << "total time: " << elapsed_time.count() << std::endl
            << std::endl;

  // Expected result
  // total time: 0.404114 (or 0.328 on new computer)

  start = std::chrono::system_clock::now();
  theta_pred = laplace_poisson_log_rng(y, n_samples, ye,
                                       sqr_exp_kernel_functor(),
                                       phi, x, delta, delta_int,
                                       theta_0, rng);
  end = std::chrono::system_clock::now();
  elapsed_time = end - start;

  std::cout << "LAPLACE_APPROX_POISSON_RNG" << std::endl
            << "total time: " << elapsed_time.count() << std::endl
            << std::endl;
}

TEST_F(laplace_disease_map_test, lk_autodiff) {
  using stan::math::var;
  using stan::math::laplace_marginal_density;
  using stan::math::diff_likelihood;

  diff_likelihood<poisson_log_likelihood> diff_functor(f, delta_lk, n_samples);

  auto start = std::chrono::system_clock::now();
  int hessian_block_size = 0;  // 0, 1, 911
  int compute_W_root = 1;
  double marginal_density_dbl
    = laplace_marginal_density(diff_functor,
                               sqr_exp_kernel_functor(),
                               value_of(phi), value_of(eta_dummy),
                               x, delta, delta_int, theta_0,
                               0, 1e-6, 100, hessian_block_size,
                               compute_W_root);

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = end - start;

  std::cout << "LAPLACE MARGINAL (dbl)" << std::endl
            << "hessian block size: " << hessian_block_size << std::endl
            << "density: " << marginal_density_dbl << std::endl
            << "total time: " << elapsed_time.count() << std::endl;

  start = std::chrono::system_clock::now();

  var marginal_density
    = laplace_marginal_density(diff_functor,
                               sqr_exp_kernel_functor(), phi, eta_dummy,
                               x, delta, delta_int, theta_0,
                               0, 1e-6, 100, hessian_block_size);

  end = std::chrono::system_clock::now();
  elapsed_time = end - start;

  VEC g;
  AVEC parm_vec = createAVEC(phi(0), phi(1));
  marginal_density.grad(parm_vec, g);

  std::cout << "LAPLACE MARGINAL AND VARI CLASS" << std::endl
            << "density: " << value_of(marginal_density) << std::endl
            << "autodiff grad: " << g[0] << " " << g[1] << std::endl
            << "total time: " << elapsed_time.count() << std::endl
            << std::endl;
  // Should return consistent evaluation of density and gradient as
  // previous iteration.
  // Expected run time: 0.39 s
}

TEST_F(laplace_disease_map_test, finite_diff_benchmark) {
  ///////////////////////////////////////////////////////////////////
  // finite_diff benchmark
  using stan::math::var;
  using stan::math::laplace_marginal_density;
  using stan::math::diff_likelihood;

  diff_likelihood<poisson_log_likelihood> diff_functor(f, delta_lk, n_samples);

  Eigen::VectorXd phi_dbl = value_of(phi);
  Eigen::VectorXd phi_u0 = phi_dbl, phi_u1 = phi_dbl,
                  phi_l0 = phi_dbl, phi_l1 = phi_dbl;
  double eps = 1e-7;

  int hessian_block_size = 1;

  phi_u0(0) += eps;
  phi_u1(1) += eps;
  phi_l0(0) -= eps;
  phi_l1(1) -= eps;

  double target_u0 = laplace_marginal_density(diff_functor,
                                 sqr_exp_kernel_functor(),
                                 phi_u0, value_of(eta_dummy),
                                 x, delta, delta_int, theta_0,
                                 0, 1e-6, 100, hessian_block_size),

  target_u1 = laplace_marginal_density(diff_functor,
                                 sqr_exp_kernel_functor(),
                                 phi_u1, value_of(eta_dummy),
                                 x, delta, delta_int, theta_0,
                                 0, 1e-6, 100, hessian_block_size),

  target_l0 = laplace_marginal_density(diff_functor,
                                       sqr_exp_kernel_functor(),
                                       phi_l0, value_of(eta_dummy),
                                       x, delta, delta_int, theta_0,
                                       0, 1e-6, 100, hessian_block_size),

  target_l1 = laplace_marginal_density(diff_functor,
                                       sqr_exp_kernel_functor(),
                                       phi_l1, value_of(eta_dummy),
                                       x, delta, delta_int, theta_0,
                                       0, 1e-6, 100, hessian_block_size);

  std::cout << "Finite_diff benchmark: " << std::endl
            << "grad: " << (target_u0 - target_l0) / (2 * eps)
            << " " << (target_u1 - target_l1) / (2 * eps)
            << std::endl;
}

TEST_F(laplace_disease_map_test, rng_autodiff) {
  using stan::math::var;
  using stan::math::laplace_base_rng;
  using stan::math::diff_likelihood;

  diff_likelihood<poisson_log_likelihood> diff_functor(f, delta_lk, n_samples);

  boost::random::mt19937 rng;
  int hessian_block_size = 0;
  int compute_W_root = 1;

  auto start = std::chrono::system_clock::now();
  Eigen::VectorXd
    theta_pred = laplace_base_rng(diff_functor,
                                  sqr_exp_kernel_functor(),
                                  phi, eta_dummy,
                                  x, x, delta, delta_int, theta_0, rng,
                                  0, 1e-6, 100, hessian_block_size,
                                  compute_W_root);
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = end - start;
  std::cout << "LAPLACE_APPROX_RNG" << std::endl
            << "total time: " << elapsed_time.count() << std::endl
            << std::endl;
}

TEST_F(laplace_disease_map_test, lpmf_wrapper) {
  using stan::math::var;
  using stan::math::laplace_marginal_lpmf;

  int hessian_block_size = 0;
  int compute_W_root = 1;

  var marginal_density
    = laplace_marginal_lpmf(n_samples, poisson_log_likelihood(),
                            eta_dummy, delta_lk,
                            sqr_exp_kernel_functor(),
                            phi, x, delta, delta_int, theta_0);

  VEC g;
  AVEC parm_vec = createAVEC(phi(0), phi(1));
  marginal_density.grad(parm_vec, g);

  std::cout << "LAPLACE MARGINAL LPMF AND VARI CLASS" << std::endl
            << "density: " << value_of(marginal_density) << std::endl
            << "autodiff grad: " << g[0] << " " << g[1] << std::endl
            << std::endl;
}

TEST_F(laplace_disease_map_test, rng_wrapper) {
  using stan::math::var;
  using stan::math::laplace_rng;

  // TODO: put these variables in the test class.
  boost::random::mt19937 rng;
  int hessian_block_size = 0;
  int compute_W_root = 1;

  Eigen::VectorXd
    theta_pred = laplace_rng(poisson_log_likelihood(),
                             eta_dummy, delta_lk, n_samples,
                             sqr_exp_kernel_functor(),
                             phi, x, delta, delta_int, theta_0, rng);

  // std::cout << "theta_pred: " << theta_pred.transpose().head(5) << std::endl;

}
