#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <stan/math/prim/mat/fun/algebra_solver.hpp>

// will want to remove these eventually
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::Dynamic;

template <typename T>
inline
Eigen::Matrix<T, Eigen::Dynamic, 1>
algebraEq(const Eigen::VectorXd x,
          const std::vector<T>& parms,
          const std::vector<double>& dat,
          const std::vector<int>& dat_int) {
  T a0 = parms[0], a1 = parms[1];
  Eigen::Matrix<T, Dynamic, 1> y(2);
  y(0) = x(0) - a0;
  y(1) = x(1) - a1;
  return y;
}

template <typename T>
struct algebraEq_functor {
private:
  std::vector<T> parms;
  std::vector<double> dat;
  std::vector<int> dat_int;
  
public:
  algebraEq_functor(const std::vector<T> parms_para,
                    const std::vector<double> dat_para,
                    const std::vector<int> dat_int_para) {
    parms = parms_para;
    dat = dat_para;
    dat_int = dat_int_para;
  }

  inline Eigen::VectorXd
  operator()(const Eigen::VectorXd x) const {
    return algebraEq(x, parms, dat, dat_int);
  }
};

template <typename T>
inline Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
jacobian(const Eigen::VectorXd x,
         const std::vector<T>& parms,
         const std::vector<double>& dat,
         const std::vector<int>& dat_int) {
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> y(2, 2);
  y(0, 0) = 1;
  y(0, 1) = 0;
  y(1, 0) = 0;
  y(1, 1) = 1;
  return y;
}

template <typename T>
struct jacobian_functor {
private:
  std::vector<T> parms;
  std::vector<double> dat;
  std::vector<int> dat_int;

public:
  jacobian_functor(const std::vector<T>& parms_para,
                   const std::vector<double>& dat_para,
                   const std::vector<int>& dat_int_para) {
    parms = parms_para;
    dat = dat_para;
    dat_int = dat_int_para;
  }

  inline
  Eigen::MatrixXd
  operator()(const Eigen::VectorXd x) const {
    return jacobian(x, parms, dat, dat_int);
  }
};


TEST(MathMatrix, dogleg_eq1) {
  Eigen::VectorXd x(2);
  x << 32, 5;

  std::vector<double> parms(2);
  parms[0] = 36;
  parms[1] = 6;

  std::vector<double> dummy_dat(0);
  std::vector<int> dummy_dat_int(0);

  algebraEq_functor<double> equation(parms, dummy_dat, dummy_dat_int);
  jacobian_functor<double> jacobian(parms, dummy_dat, dummy_dat_int);

  Eigen::VectorXd theta;
  theta = stan::math::algebra_solver(equation, 
                                     jacobian,
                                     x, parms, dummy_dat, dummy_dat_int);

  EXPECT_EQ(36, theta(0));
  EXPECT_EQ(6, theta(1));
}
