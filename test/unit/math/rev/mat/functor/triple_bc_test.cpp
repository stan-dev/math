#include <stan/math/rev/mat.hpp>
#include <stan/math/rev/mat/functor/triple.hpp>
#include <gtest/gtest.h>
#include <iostream>

TEST(MathMatrix, triple3) {
  using stan::math::var;
  var a = 10;
  var b = 20;

  std::vector<double> g(3);
  std::vector<var> ab(2);
  ab[0] = a;
  ab[1] = b;
  Eigen::Matrix<var, Eigen::Dynamic, 1> x(2);
  x(0) = a;
  x(1) = b;

  Eigen::Matrix<var, Eigen::Dynamic, 1> c = stan::math::triple(x);

  std::cout << "This is the last thing you'll see before the code breaks..."
            << std::endl;
  c(0).grad(ab, g);
  std::cout << "c(0) = " << c(0) << "; g.size() = " << g.size()
            << "; g[0] = " << g[0] << "; g[1] = " << g[1] << std::endl;
}
