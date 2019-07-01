
#include <stan/math/rev/core.hpp>
#include <stan/math/laplace/lgp_solver.hpp>
#include <stan/math/rev/mat/functor/algebra_system.hpp>
#include <stan/math/rev/mat/functor/kinsol_data.hpp>
#include <stan/math/rev/mat/functor/kinsol_solve.hpp>

#include <test/unit/math/rev/mat/functor/util_algebra_solver.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>

TEST(laplace, basics) {
  // Simple experiment using Eigen's map.
  // Map creates a pointer, but a matirx set equal to the
  // map will cause a duplicate of the data.
  
  int N = 3;
  int dim = 2 * N + N * N;

  // double n_samples_arr[dim];
  std::vector<double> n_samples_arr(dim);
  for (int i = 0; i < dim; i++) n_samples_arr[i] = i;

  Eigen::Map<Eigen::VectorXd> n_samples(&n_samples_arr[0], N);
  Eigen::Map<Eigen::VectorXd> sums(&n_samples_arr[N], N);
  Eigen::Map<Eigen::MatrixXd> Q(&n_samples_arr[2 * N], N, N);
  
  Eigen::MatrixXd Q_matrix = Q;

  // std::cout << n_samples.transpose() << std::endl;
  // std::cout << sums.transpose() << std::endl << std::endl;
  // std::cout << Q_matrix << std::endl << std::endl;

  // Test to show Q_matrix is a duplicate (though Q is not)
  n_samples_arr[14] = -1;
  // std::cout << Q_matrix << std::endl;
}

TEST(laplace, system) {
  using stan::math::lgp_dense_system;
  using stan::math::to_array_1d;
  using stan::math::lgp_f;
  using stan::math::lgp_J_f;
  
  int dim_theta = 2;
  Eigen::VectorXd theta(2);
  theta << -0.7268203, 1.3347728;
  
  Eigen::VectorXd phi(2);
  phi << 4, 0;
  
  Eigen::VectorXd n_samples(dim_theta);
  for (int i = 0; i < dim_theta; i++) n_samples(i) = 5;
  Eigen::VectorXd sums(dim_theta);
  sums << 3, 10;

  // TO DO -- efficient way to construct the data vector.
  // (actually, should be straightfoward)
  lgp_dense_system<double> system(phi, n_samples, sums);

  int n_elements = dim_theta * (2 + dim_theta);
  std::vector<double> dat(n_elements);
  for (int i = 0; i < dim_theta; i++) dat[i] = n_samples[i];
  for (int i = 0; i < dim_theta; i++) dat[dim_theta + i] = sums[i];

  std::vector<double> Q_array = to_array_1d(system.Q_);
  for (int i = 0; i < dim_theta * dim_theta; i++)
    dat[2 * dim_theta + i] = Q_array[i];
  std::vector<int> dat_int;
  
  lgp_f gradient_eval;
  Eigen::VectorXd f = gradient_eval(theta, phi, dat, dat_int, 0);

  EXPECT_FLOAT_EQ(0.7644863, f(0));
  EXPECT_FLOAT_EQ(-9.3293567, f(1));

  // std::cout << gradient_eval(theta, phi, dat, dat_int, 0) << std::endl << std::endl;

  lgp_J_f hessian_eval;
  double x_sun[2];
  x_sun[0] = theta(0);
  x_sun[1] = theta(1);

  SUNMatrix J = SUNDenseMatrix(dim_theta, dim_theta);

  hessian_eval(0, theta, phi, dat, dat_int, 0, x_sun, J);

  EXPECT_FLOAT_EQ(-2.667219, SM_ELEMENT_D(J, 0, 0));
  EXPECT_FLOAT_EQ(0, SM_ELEMENT_D(J, 1, 0));
  EXPECT_FLOAT_EQ(0, SM_ELEMENT_D(J, 0, 1));
  EXPECT_FLOAT_EQ(-19.245664, SM_ELEMENT_D(J, 1, 1));
}
