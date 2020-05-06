#include <iostream>
#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathRev, TestVarEigen) {
  using stan::math::var_type;
  using stan::math::var;
  using stan::math::sum;
  Eigen::Matrix<double, -1, -1> x_vals(3, 3);
  Eigen::Matrix<double, -1, -1> y_vals(3, 3);
  x_vals << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  y_vals << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  var_type<Eigen::Matrix<double, -1, -1>> x = x_vals;
  var_type<Eigen::Matrix<double, -1, -1>> y = y_vals;
  var lp = 0;
  auto mul_xy = x * y;
  auto sum_mul_xy = sum(mul_xy + x);
  lp -= sum_mul_xy;
  lp.grad();
  puts("-------------");
  std::cout << "lp static val: \n" << lp.val() << "\n";
  std::cout << "lp static adj: \n" << lp.adj() << "\n";
  puts("---------");
  std::cout << "sum_mul_xy static val: \n" << sum_mul_xy.val() << "\n";
  std::cout << "sum_mul_xy static adj: \n" << sum_mul_xy.adj() << "\n";
  puts("---------");

  std::cout << "mul_xy stat val: \n" << mul_xy.val() << "\n";
  std::cout << "mul_xy stat adj: \n" << mul_xy.adj() << "\n";
  puts("---------");

  std::cout << "x stat val: \n" << x.val() << "\n";
  std::cout << "x stat adj: \n" << x.adj() << "\n";
  puts("---------");

  std::cout << "y stat val: \n" << y.val() << "\n";
  std::cout << "y stat adj: \n" << y.adj() << "\n";
  stan::math::set_zero_all_adjoints();
  puts("---------");
  puts("Dynamic Matrix");
  puts("---------");
  Eigen::Matrix<var, -1, -1> x_dyn = x_vals;
  Eigen::Matrix<var, -1, -1> y_dyn = y_vals;
  Eigen::Matrix<var, -1, -1> mul_xy_dyn = x_dyn * y_dyn;
  var sum_mul_xy_dyn = sum(mul_xy_dyn + x_dyn);
  var lp_dyn = 0;
  lp_dyn -= sum_mul_xy_dyn;
  lp.grad();
  std::cout << "lp dyn val: \n" << lp.val() << "\n";
  std::cout << "lp dyn adj: \n" << lp.adj() << "\n";
  puts("---------");
  std::cout << "sum_mul_xy dyn val: \n" << sum_mul_xy.val() << "\n";
  std::cout << "sum_mul_xy dyn adj: \n" << sum_mul_xy.adj() << "\n";
  puts("---------");

  std::cout << "mul_xy dyn val: \n" << mul_xy.val() << "\n";
  std::cout << "mul_xy dyn adj: \n" << mul_xy.adj() << "\n";
  puts("---------");

  std::cout << "x dyn val: \n" << x.val() << "\n";
  std::cout << "x dyn adj: \n" << x.adj() << "\n";
  puts("---------");

  std::cout << "y dyn val: \n" << y.val() << "\n";
  std::cout << "y dyn adj: \n" << y.adj() << "\n";
}
