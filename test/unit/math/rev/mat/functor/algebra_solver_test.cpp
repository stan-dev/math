#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/functor/algebra_solver.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <gtest/gtest.h>
#include <iostream>


template <typename T1, typename T2>
inline
Eigen::Matrix<typename boost::math::tools::promote_args<T1, T2>::type, 
              Eigen::Dynamic, 1>
algebraEq(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
          const Eigen::Matrix<T2, Eigen::Dynamic, 1>& parms,
          const std::vector<double>& dat,
          const std::vector<int>& dat_int) {
  typedef typename boost::math::tools::promote_args<T1, T2>::type scalar;
  Eigen::Matrix<scalar, Eigen::Dynamic, 1> y(2);
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
  using stan::math::algebra_solver;
  using stan::math::sum;
  using stan::math::var;

  for (size_t k = 0; k < 2; k++) {
    Eigen::VectorXd x(2);
    x << 1, 1;  // initial guess
    
    Eigen::Matrix<var, Eigen::Dynamic, 1> parms(2);
    parms << 36, 6;

    std::vector<double> dummy_dat(0);
    std::vector<int> dummy_dat_int(0);

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta;
    theta = algebra_solver(algebraEq_functor<double, double>(),
                           x, parms, dummy_dat, dummy_dat_int);
    EXPECT_EQ(36, theta(0));
    EXPECT_EQ(6, theta(1));

    std::vector<double> dtheta_dparms(parms.size());
    if (k == 0) {
      dtheta_dparms[0] = 1;
      dtheta_dparms[1] = 0;
    }
    else {
      dtheta_dparms[0] = 0;
      dtheta_dparms[1] = 1;
    }

    AVEC parms_vec = createAVEC(parms(0), parms(1));
    VEC g;
    theta(k).grad(parms_vec, g);

    for (int i = 0; i < parms.size(); i++)
      EXPECT_EQ(dtheta_dparms[i], g[i]);
  }
  
}

/*
template <typename T1, typename T2>
inline Eigen::Matrix<typename boost::math::tools::promote_args<T1, T2>::type,
                     Eigen::Dynamic, 1>
algebraEq2(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
           const Eigen::Matrix<T2, Eigen::Dynamic, 1>& parms,
           const std::vector<double>& dat,
           const std::vector<int>& dat_int) {
  typedef typename boost::math::tools::promote_args<T1, T2>::type scalar;
  Eigen::Matrix<scalar, Eigen::Dynamic, 1> y(3);
  y(0) = x(0) - parms(0);
  y(1) = x(1) - parms(1);
  y(2) = x(2) - parms(2);
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
      return algebraEq(x, parms, dat, dat_int);
    else
      return algebraEq(theta, x, dat, dat_int);
  }
};

TEST(MathMatrix, algebra_solver_eq2) {
  using stan::math::algebra_solver;
  using stan::math::sum;
  using stan::math::var;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  
  int nParms = 3, nSolutions = 3;
  for (int k = 0; k < nParms; k++) {
    Eigen::VectorXd x(nSolutions);
    x << 1, 1, 1;  // initial guess
    
    Matrix<var, Dynamic, 1> parms(nParms);
    parms << 36, 6, 42;
    
    std::vector<double> dummy_dat(0);
    std::vector<int> dummy_dat_int(0);
    
    Matrix<var, Dynamic, 1> theta;
    theta = algebra_solver(algebraEq2_functor<double, double>(),
                           x, parms, dummy_dat, dummy_dat_int);
    
    EXPECT_EQ(36, theta(0));
    EXPECT_EQ(6, theta(1));
    EXPECT_EQ(42, theta(2));
    
    Matrix<double, Dynamic, Dynamic> dtheta_dparms(nParms, nSolutions);
    dtheta_dparms << 1, 0, 0,
                     0, 1, 0,
                     0, 0, 1;
    
    AVEC parms_vec = createAVEC(parms(0), parms(1), parms(2));
    VEC g;
    theta(k).grad(parms_vec, g);
    for (int i = 0; i < nParms; i++)
      EXPECT_EQ(dtheta_dparms(k, i), g[i]);
  }
} */
