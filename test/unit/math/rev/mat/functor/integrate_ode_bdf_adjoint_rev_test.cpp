#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <boost/numeric/odeint.hpp>
#include <test/unit/math/rev/mat/functor/util_cvodes_bdf_adjoint.hpp>
#include <test/unit/math/prim/arr/functor/harmonic_oscillator.hpp>
#include <test/unit/math/prim/arr/functor/lorenz.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <chrono>

TEST(StanAgradRevOde_integrate_ode, harmonic_oscillator_basic) {
  stan::math::ChainableStack::instance().var_stack_.reserve(10000);

  using stan::math::var;
  harm_osc_ode_fun harm_osc;

  /*std::vector<double> theta;
  theta.push_back(0.15);

  std::vector<double> y0;
  y0.push_back(1.0);
  y0.push_back(0.0);*/
  lorenz_ode_fun lorenz;

  std::vector<double> y0;
  std::vector<double> theta;
  double t0;

  t0 = 0;

  theta.push_back(10.0);
  theta.push_back(28.0);
  theta.push_back(8.0 / 3.0);
  y0.push_back(10.0);
  y0.push_back(1.0);
  y0.push_back(1.0);

  std::vector<double> ts;
  for (int i = 0; i < 100; i++)
    ts.push_back(t0 + 0.1 * (i + 1));

  std::vector<double> x;
  std::vector<int> x_int;

  using stan::math::promote_scalar;
  using stan::math::var;

  {
    auto y0_ = promote_scalar<var>(y0);
    auto theta_ = promote_scalar<var>(theta);
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<std::vector<var> > ode_res_vd
        = stan::math::integrate_ode_bdf_adjoint(lorenz, y0_, t0, ts, theta_, x,
                                                x_int);
    auto middle = std::chrono::high_resolution_clock::now();

    var z = 0.0;
    for (int i = 0; i < ode_res_vd.size(); ++i)
      for (int j = 0; j < ode_res_vd[i].size(); ++j)
        z += ode_res_vd[i][j];

    z.grad();

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "adjoint" << std::endl;
    std::cout << "forward : "
              << std::chrono::duration<double>(middle - start).count()
              << std::endl;
    std::cout << "reverse : "
              << std::chrono::duration<double>(end - middle).count()
              << std::endl;
    std::cout << "time : " << std::chrono::duration<double>(end - start).count()
              << std::endl;
    std::cout << z.val() << std::endl;
    for (int i = 0; i < y0_.size(); ++i)
      std::cout << y0_[i].adj() << std::endl;
    for (int i = 0; i < theta_.size(); ++i)
      std::cout << theta_[i].adj() << std::endl;
  }

  stan::math::recover_memory();

  {
    auto y0_ = promote_scalar<var>(y0);
    auto theta_ = promote_scalar<var>(theta);
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<std::vector<var> > ode_res_vd
        = stan::math::integrate_ode_bdf(lorenz, y0_, t0, ts, theta_, x, x_int);

    var z = 0.0;
    for (int i = 0; i < ode_res_vd.size(); ++i)
      for (int j = 0; j < ode_res_vd[i].size(); ++j)
        z += ode_res_vd[i][j];

    z.grad();
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "regular" << std::endl;
    std::cout << "time : " << std::chrono::duration<double>(end - start).count()
              << std::endl;
    std::cout << z.val() << std::endl;
    for (int i = 0; i < y0_.size(); ++i)
      std::cout << y0_[i].adj() << std::endl;
    for (int i = 0; i < theta_.size(); ++i)
      std::cout << theta_[i].adj() << std::endl;
  }

  stan::math::recover_memory();
}
