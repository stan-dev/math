#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathRev, TestVarEigen) {
  using stan::math::var_type;
  using stan::math::var;
  using stan::math::sum;
  Eigen::Matrix<double, -1, -1> x_vals = Eigen::MatrixXd::Random(5, 5);
  Eigen::Matrix<double, -1, -1> y_vals = Eigen::MatrixXd::Random(5, 5);
  var_type<Eigen::Matrix<double, -1, -1>> x = x_vals;
  var_type<Eigen::Matrix<double, -1, -1>> y = y_vals;
  var_type<double> lp = 0;
  lp -= sum(x * y);
  puts("Static Matrix:");
  std::cout << "static vals: \n" << lp.val() << "\n";
  lp.grad();
  std::cout << "static x adj: \n" << x.adj() << "\n";
  std::cout << "static y adj: \n" << y.adj() << "\n";
  puts("---------");


  puts("Dynamic Matrix:");
  Eigen::Matrix<var, -1, -1> x_dyn = x_vals;
  Eigen::Matrix<var, -1, -1> y_dyn = y_vals;
  var_type<double> lp_dyn = 0;
  lp_dyn -= sum(x_dyn * y_dyn);
  std::cout << "dynamic vals: \n" << lp_dyn.val() << "\n";
  lp_dyn.grad();
  std::cout << "dynamic x adj: \n" << x_dyn.adj() << "\n";
  std::cout << "dynamic y adj: \n" << y_dyn.adj() << "\n";
/*
  using std::pow;
  double yy = 1.3;
  stan::math::var mu = 0.5, sigma = 1.2;
  stan::math::var lp_n = 0;
  lp_n -= 0.5 * log(2 * stan::math::pi());
  lp_n -= log(sigma);
  lp_n -= 0.5 * pow((yy - mu) / sigma, 2);
  std::cout << "f(mu, sigma) = " << lp_n.val() << std::endl;
  lp_n.grad();
  std::cout << " d.f / d.mu = " << mu.adj()
   << " d.f / d.sigma = " << sigma.adj() << std::endl;
*/}
