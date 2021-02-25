#include <stan/math.hpp>
#include <stan/math/laplace/laplace.hpp>
#include <stan/math/laplace/laplace_likelihood_general.hpp>
// #include <stan/math/laplace/laplace_likelihood_poisson_log.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>

#include <test/unit/math/laplace/laplace_utility.hpp>
#include <test/unit/math/rev/fun/util.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>

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
};


TEST_F(laplace_disease_map_test, lk_analytical) {
  // Based on (Vanhatalo, Pietilainen and Vethari, 2010). See
  // https://research.cs.aalto.fi/pml/software/gpstuff/demo_spatial1.shtml
  using stan::math::var;
  using stan::math::laplace_marginal_poisson_log_lpmf;
  using stan::math::sqr_exp_kernel_functor;

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
  /*
  using stan::math::diff_poisson_log;
  using stan::math::to_vector;
  using stan::math::sqr_exp_kernel_functor;
  using stan::math::laplace_rng;
  using stan::math::laplace_poisson_log_rng;

  diff_poisson_log diff_likelihood(to_vector(n_samples),
                                   to_vector(y),
                                   stan::math::log(ye));
  boost::random::mt19937 rng;
  start = std::chrono::system_clock::now();
  Eigen::VectorXd
    theta_pred = laplace_rng(diff_likelihood,
                            sqr_exp_kernel_functor(),
                            phi, x, delta, delta_int,
                            theta_0, rng);

  end = std::chrono::system_clock::now();
  elapsed_time = end - start;

  std::cout << "LAPLACE_APPROX_RNG" << std::endl
            << "total time: " << elapsed_time.count() << std::endl
            << std::endl;

  // Expected result
  // total time: 0.404114

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
  */
}

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

TEST_F(laplace_disease_map_test, lk_autodiff) {
  using stan::math::var;
  using stan::math::laplace_marginal_density;
  using stan::math::diff_likelihood;
  using stan::math::sqr_exp_kernel_functor;

  Eigen::VectorXd delta_lk(2 * n_observations);
  for (int i = 0; i < n_observations; i++) delta_lk(i) = y[i];
  for (int i = 0; i < n_observations; i++) delta_lk(n_observations + i) = ye(i);

  poisson_log_likelihood f;
  diff_likelihood<poisson_log_likelihood>
    diff_functor(f, delta_lk, n_samples);

  auto start = std::chrono::system_clock::now();

  Eigen::Matrix<var, -1, 1> eta_dummy;
  var marginal_density
    = laplace_marginal_density(diff_functor,
                               sqr_exp_kernel_functor(), phi, eta_dummy,
                               x, delta, delta_int, theta_0);

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
  // Should return consistent evaluation of density and gradient as
  // previous iteration.
  // Expected run time: 0.39 s
}
