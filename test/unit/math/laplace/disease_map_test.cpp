#include <stan/math.hpp>
#include <stan/math/laplace/laplace_marginal_poisson.hpp>

#include <test/unit/math/laplace/lgp_utility.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>


TEST(laplace, disease_map_dim_911) {
  // Based on (Vanhatalo, Pietilainen and Vethari, 2010). See
  // https://research.cs.aalto.fi/pml/software/gpstuff/demo_spatial1.shtml
  using stan::math::var;
  using stan::math::laplace_marginal_poisson;

  int dim_theta = 911;
  int n_observations = 911;
  std::string data_directory = "test/unit/math/laplace/aki_disease_data/";
  std::vector<double> x1(dim_theta), x2(dim_theta);
  std::vector<int> y(n_observations);
  Eigen::VectorXd ye(n_observations);
  read_in_data(dim_theta, n_observations, data_directory, x1, x2, y, ye);

  // look at some of the data
  std::cout << "x_1: " << x1[0] << " " << x2[0] << std::endl
            << "x_2: " << x1[1] << " " << x2[1] << std::endl
            << "y_1: " << y[0] << " y_2: " << y[1] << std::endl
            << "ye_1: " << ye[0] << " ye_2: " << ye[1] << std::endl;

  int dim_x = 2;
  std::vector<Eigen::VectorXd> x(dim_theta);
  for (int i = 0; i < dim_theta; i++) {
    Eigen::VectorXd coordinate(dim_x);
    coordinate << x1[i], x2[i];
    x[i] = coordinate;
  }

  // one observation per group
  std::vector<int> n_samples(dim_theta);
  for (int i = 0; i < dim_theta; i++) n_samples[i] = 1;

  Eigen::VectorXd theta_0 = Eigen::VectorXd::Zero(dim_theta);
  int dim_phi = 2;
  Eigen::Matrix<var, Eigen::Dynamic, 1> phi(dim_phi);
  phi << 0.3162278, 200;  // variance, length scale

  auto start = std::chrono::system_clock::now();

  var marginal_density = laplace_marginal_poisson(theta_0, phi, x,
                                                  n_samples, y, ye);

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
}
