#include <gtest/gtest.h>

#include <iostream>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <ctime>

#include <boost/random/mersenne_twister.hpp>

#include <stan/math.hpp>

// very small michaelis menten example
#include <test/unit/math/rev/arr/functor/coupled_mm.hpp>


TEST(StanOde_stiff_small_stress_test, cvode_coupled_mm) {

  const std::clock_t clock_start = std::clock();

  const size_t nIntegrations = 100;
  const double rel_tol = 1E-10;
  const double abs_tol = 1E-10;
  
  coupled_mm_ode_fun f_;

  // initial value and parameters from model definition
  std::vector<double> y0(2);
  y0[0] = 1.0;
  y0[1] = 1E-3;

  double t0 = 0;

  std::vector<double> ts;
  const double tmax = 100;

  for (int i = 0; i < 101; i++)
    ts.push_back((i+1)*tmax/101.);

  std::vector<double> theta(4);

  theta[0] = 1.0;
  theta[1] = 0.5;
  theta[2] = 0.5;
  theta[3] = 0.1;

  boost::random::mt19937 rng;

  std::vector<double> data;

  std::vector<int> data_int;

  size_t integration = 0;

  rng.seed(45656);
  for( ; integration < nIntegrations; integration++) {

    std::vector<double> theta_run(theta);

    for(size_t i = 0; i < theta.size(); i++)
      theta_run[i] *= stan::math::lognormal_rng(0, 3.0, rng);

    std::vector<double> y0_run(y0);

    for(size_t i = 0; i < y0.size(); i++)
      y0_run[i] *= stan::math::lognormal_rng(0, 3.0, rng);

    std::vector<stan::math::var> theta_var(theta_run.begin(), theta_run.end());
    std::vector<stan::math::var> y0_var(y0_run.begin(), y0_run.end());

    std::vector<std::vector<stan::math::var> > res_cvodes
     = stan::math::integrate_ode_cvodes(f_, y0_var, t0, ts, theta_var, data, data_int, rel_tol, abs_tol, 1E10, 1);
    stan::math::recover_memory();

    std::clock_t clock_end = std::clock();
    double run_duration_sec = (clock_end - clock_start) / CLOCKS_PER_SEC ;

    if(run_duration_sec > 60)
      break;

  }

  EXPECT_EQ(integration, nIntegrations);
 }
