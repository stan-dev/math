#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/propogate_static_matrix.hpp>
#include <gtest/gtest.h>

TEST(RevCore, static_to_dynamic_multi_grad) {
  using stan::math::var_value;
  Eigen::MatrixXd v = Eigen::MatrixXd::Ones(3, 3);

  var_value<Eigen::MatrixXd> x1(v);
  Eigen::Matrix<var_value<double>, -1, -1> x2 = x1;
  puts("Grad x2(0, 0): ");
  x2(0, 0).grad();
  std::cout << x1.adj() << std::endl;
  puts("Grad x2(1, 0): ");
  x2(1, 0).grad();
  std::cout << x1.adj() << std::endl;
  puts("Grad x2(1, 1): ");
  x2(1, 1).grad();
  std::cout << x1.adj() << std::endl;
  x2(0, 1) = var_value<double>(10);
  x2(0, 1).grad();
  puts("Grad x2(0, 1) new: ");
  std::cout << x1.adj() << std::endl;

}

TEST(RevCore, static_to_dynamic_multi_grad_clear) {
  using stan::math::var_value;
  Eigen::MatrixXd v = Eigen::MatrixXd::Ones(3, 3);

  var_value<Eigen::MatrixXd> x1(v);
  Eigen::Matrix<var_value<double>, -1, -1> x2 = x1;
  puts("Grad x2(0, 0): ");
  x2(0, 0).grad();
  std::cout << x1.adj() << std::endl;
  stan::math::set_zero_all_adjoints();
  puts("Grad x2(1, 0): ");
  x2(1, 0).grad();
  std::cout << x1.adj() << std::endl;
  stan::math::set_zero_all_adjoints();
  puts("Grad x2(1, 1): ");
  x2(1, 1).grad();
  std::cout << x1.adj() << std::endl;
  stan::math::set_zero_all_adjoints();
}

TEST(RevCore, dynamic_to_static) {
  using stan::math::var_value;
  Eigen::MatrixXd v = Eigen::MatrixXd::Ones(3, 3);
  Eigen::Matrix<var_value<double>, -1, -1> x1(v);
  var_value<Eigen::MatrixXd> x2(x1);
  x1(1, 1) = var_value<double>(10);
  x2.grad();
  puts("X2: ");
  std::cout << x2.adj() << std::endl;
  puts("X1: ");
  std::cout << x1.adj() << std::endl;
}

TEST(RevCore, dynamic_to_dynamic) {
  using stan::math::var_value;
  Eigen::MatrixXd v = Eigen::MatrixXd::Ones(3, 3);
  Eigen::Matrix<var_value<double>, -1, -1> x1(v);
  Eigen::Matrix<var_value<double>, -1, -1> x2 = x1;
  x1(1, 1) = var_value<double>(10);
  for (int i = 0; i < x2.size(); ++i) {
    x2(i).grad();
  }
  puts("X2: ");
  std::cout << x2.adj() << std::endl;
  puts("X1: ");
  std::cout << x1.adj() << std::endl;
}
