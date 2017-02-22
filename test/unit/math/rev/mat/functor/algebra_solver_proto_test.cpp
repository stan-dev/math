#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat.hpp>  // probably don't need this include
#include <gtest/gtest.h>
#include <stan/math/rev/mat/functor/algebra_solver.hpp>

// will want to remove these eventually
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::Dynamic;

template <typename T1, typename T2>
inline
Eigen::Matrix<typename boost::math::tools::promote_args<T1, T2>::type, 
              Eigen::Dynamic, 1>
algebraEq(const Eigen::Matrix<T1, Eigen::Dynamic, 1> x,
          const Eigen::Matrix<T2, Eigen::Dynamic, 1>& parms,
          const std::vector<double>& dat,
          const std::vector<int>& dat_int) {
  typedef typename boost::math::tools::promote_args<T1, T2>::type scalar;
  Eigen::Matrix<scalar, Dynamic, 1> y(2);
  y(0) = x(0) - parms(0);
  y(1) = x(1) - parms(1);
  return y;
}

template <typename T1, typename T2>
struct algebraEq_functor {
private:
  Eigen::Matrix<T1, Eigen::Dynamic, 1> theta;
  Eigen::Matrix<T2, Eigen::Dynamic, 1> parms;
  std::vector<double> dat;
  std::vector<int> dat_int;
  std::string variable;

public:
  algebraEq_functor() { };

  algebraEq_functor(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta_para,
                    const Eigen::Matrix<T2, Eigen::Dynamic, 1>& parms_para,
                    const std::vector<double>& dat_para,
                    const std::vector<int>& dat_int_para,
                    const std::string& variable_para) {
    theta = theta_para;
    parms = parms_para;
    dat = dat_para;
    dat_int = dat_int_para;
    variable = variable_para;
  }

  template <typename T>
  inline
  Eigen::Matrix<typename boost::math::tools::promote_args<T, T1, T2>::type,
                Eigen::Dynamic, 1>
  operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1> x) const {
    if (variable == "theta")
      return algebraEq(x, parms, dat, dat_int);
    else
      return algebraEq(theta, x, dat, dat_int);
  }
};


TEST(MathMatrix, algebra_solver_eq1) {
  Eigen::VectorXd x(2), parms(2);
  x << 32, 5;
  parms << 36, 6;

  std::vector<double> dummy_dat(0);
  std::vector<int> dummy_dat_int(0);

  Eigen::VectorXd theta;
  theta = stan::math::algebra_solver(algebraEq_functor<double, double>(),
                                     x, parms, dummy_dat, dummy_dat_int);

  EXPECT_EQ(36, theta(0));
  EXPECT_EQ(6, theta(1));
}


template <typename T1, typename T2>
inline
Eigen::Matrix<typename boost::math::tools::promote_args<T1, T2>::type, 
              Eigen::Dynamic, 1>
algebraEq2(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
           const Eigen::Matrix<T2, Eigen::Dynamic, 1>& parms,
           const std::vector<double>& dat,
           const std::vector<int>& dat_int) {
  typedef typename boost::math::tools::promote_args<T1, T2>::type scalar;
  Eigen::Matrix<scalar, Eigen::Dynamic, 1> y(3);
  y(0) = x(0) - parms(0);
  y(1) = x(1) - parms(1);
  y(2) = x(2) * parms(2);
  return y;
}

template <typename T1, typename T2>
struct algebraEq2_functor {
private:
  Eigen::Matrix<T1, Eigen::Dynamic, 1> theta;
  Eigen::Matrix<T2, Eigen::Dynamic, 1> parms;
  std::vector<double> dat;
  std::vector<int> dat_int;
  std::string variable;

public:
  algebraEq2_functor() { };

  algebraEq2_functor(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta_para,
                     const Eigen::Matrix<T2, Eigen::Dynamic, 1>& parms_para,
                     const std::vector<double>& dat_para,
                     const std::vector<int>& dat_int_para,
                     const std::string& variable_para) {
    theta = theta_para;
    parms = parms_para;
    dat = dat_para;
    dat_int = dat_int_para;
    variable = variable_para;
  }

  template <typename T>
  inline
  Eigen::Matrix<typename boost::math::tools::promote_args<T, T1, T2>::type,
                Eigen::Dynamic, 1>
  operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1> x) const {
    if (variable == "theta")
      return algebraEq2(x, parms, dat, dat_int);
    else
      return algebraEq2(theta, x, dat, dat_int);
  }
};

TEST(MathMatrix, algebra_solver_eq2) {
  Eigen::VectorXd x(3), parms(3);
  x << 32, 5, 2;
  parms << 36, 6, 42;

  std::vector<double> dummy_dat(0);
  std::vector<int> dummy_dat_int(0);

  Eigen::VectorXd theta;
  theta = stan::math::algebra_solver(algebraEq2_functor<double, double>(),
                                     x, parms, dummy_dat, dummy_dat_int);

  EXPECT_EQ(36, theta(0));
  EXPECT_EQ(6, theta(1));
  EXPECT_EQ(0, theta(2));
}

