#include <stan/math/prim/mat.hpp>
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


// will want to remove these eventually
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::Dynamic;

/***
 * Implement original Eigen hybrid solver.
 */
template <typename _Scalar, int NX = Dynamic, int NY = Dynamic>
struct Functor {
  typedef _Scalar Scalar;
  enum {
    InputsAtCompileTime = NX,
    ValuesAtCompileTime = NY
  };
  
  typedef Matrix<Scalar, InputsAtCompileTime, 1> InputType;
  typedef Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
  typedef Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;
  
  const int m_inputs, m_values;
  
  Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
  Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}
  
  int inputs() const { return m_inputs; }
  int values() const { return m_values; }
};

struct hybrj_functor : Functor<double> {
  hybrj_functor(void) : Functor<double>(9, 9) { }
  
  int operator()(const VectorXd& x, VectorXd& fvec) {
    fvec = algebraEq(x);
    return 0;
  }

  int df(const VectorXd& x, MatrixXd& fjac) {
    fjac = jacobian(x);
    return 0;
  }
};

TEST(MathMatrix, hybrid_eigen) {

  const int n = 2;
  VectorXd x(n);
  x << 32, 5;

  // do the computation
  hybrj_functor functor;
  Eigen::HybridNonLinearSolver<hybrj_functor> solver(functor);
  Eigen::VectorXd theta = x, fvec;
  // solver(x, fvec);
  // std::cout << fvec << std::endl;
  solver.solve(theta);
  
  EXPECT_EQ(36, theta(0));
  EXPECT_EQ(6, theta(1));
}

///////////////////////////////////////////////////////////////////////////////

/*
 * Test dogleg function.
 */
TEST(MathMatrix, dogleg_eq1) {
  Eigen::VectorXd x(2);
  x << 32, 5;

  Eigen::VectorXd theta;
  theta = stan::math::dogleg(x, algebraEq_functor(), jacobian_functor());

  EXPECT_EQ(36, theta(0));
  EXPECT_EQ(6, theta(1));
}


inline Eigen::VectorXd
algebraEq2(const Eigen::VectorXd x) {
  Eigen::VectorXd y(3);
  y(0) = x(0) - 36;
  y(1) = x(1) - 6;
  y(2) = x(2) * 42;
  return y;
}

struct algebraEq_functor2 {
  inline Eigen::VectorXd
  operator()(const Eigen::VectorXd x) const {
    return algebraEq2(x);
  }
};

inline Eigen::MatrixXd
jacobian2(const Eigen::VectorXd x) {
  Eigen::MatrixXd y(3, 3);
  y << 1, 0, 0,
       0, 1, 0,
       0, 0, 1;

  return y;
}

struct jacobian_functor2 {
  inline Eigen::MatrixXd
  operator()(const Eigen::VectorXd x) const {
    return jacobian2(x);
  }
};


TEST(MathMatrix, dogleg_eq2) {
  Eigen::VectorXd x(3);
  x << 32, 5, 10;

  Eigen::VectorXd theta;
  theta = stan::math::dogleg(x, algebraEq_functor2(), jacobian_functor2());

  EXPECT_EQ(36, theta(0));
  EXPECT_EQ(6, theta(1));
  EXPECT_NEAR(0, theta(2), 1e-30);  // obtained result is not exactly 0
}
