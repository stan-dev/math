#include <stan/math/prim/mat.hpp>
#include <stan/math/prim/mat/fun/dogleg.hpp>
#include <gtest/gtest.h>


inline Eigen::VectorXd
algebraEq(const Eigen::VectorXd x) {
  Eigen::VectorXd y(2);
  y(0) = x(0) - 36;
  y(1) = x(1) - 6;
  return y;
}

struct algebraEq_functor {
  inline Eigen::VectorXd
  operator()(const Eigen::VectorXd x) const {
    return algebraEq(x);
  }
};

inline Eigen::MatrixXd
jacobian(const Eigen::VectorXd x) {
  Eigen::MatrixXd y(2, 2);
  y(0, 0) = 1;
  y(0, 1) = 0;
  y(1, 0) = 0;
  y(1, 1) = 1;
  return y;
}

struct jacobian_functor {
  inline Eigen::MatrixXd
  operator()(const Eigen::VectorXd x) const {
    return jacobian(x);
  }
};

TEST(MathMatrix, dogleg) {
  Eigen::VectorXd x(2);
  x << 32, 5;

  Eigen::VectorXd theta;
  theta = stan::math::dogleg(x, algebraEq_functor(), jacobian_functor());
}
