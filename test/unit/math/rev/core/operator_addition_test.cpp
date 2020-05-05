#include <stan/math/rev.hpp>
#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathRev, TestVarEigen) {
  using stan::math::var_type;
  using stan::math::var;
  Eigen::Matrix<double, -1, -1> x_vals = Eigen::MatrixXd::Random(5, 5);
  Eigen::Matrix<double, -1, -1> y_vals = Eigen::MatrixXd::Random(5, 5);
  var_type<Eigen::Matrix<double, -1, -1>> x(x_vals);
  var_type<Eigen::Matrix<double, -1, -1>> y(y_vals);
  auto z = x + y;
  auto zz = stan::math::sum(z);
  zz.grad();
  std::cout << "static grad vals: \n" << zz.val() << "\n";
  std::cout << "static grad adj: \n" << zz.adj() << "\n";
  Eigen::Matrix<var, -1, -1> x_dyn = x_vals;
  Eigen::Matrix<var, -1, -1> y_dyn = y_vals;

  Eigen::Matrix<var, -1, -1> z_dyn = x_dyn + y_dyn;
  var zz_dyn = z_dyn.sum();
  std::cout << "dynamic vals: \n" << zz_dyn.val() << "\n";
  std::cout << "dynamic adj: \n" << zz_dyn.adj() << "\n";
  zz_dyn.grad();
  std::cout << "dynamic grad vals: \n" << zz_dyn.val() << "\n";
  std::cout << "dynamic grad adj: \n" << zz_dyn.adj() << "\n";

}
