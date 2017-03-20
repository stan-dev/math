#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/functor/algebra_solver.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <stan/math/rev/core.hpp>
// #include <stan/math/rev/core/chainable_alloc.hpp>
// #include <stan/math/rev/core/chainablestack.hpp>
#include <gtest/gtest.h>
#include <iostream>

/* template <typename T1, typename T2>
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


TEST(MathMatrix, algebra_solver_chain) {
  using stan::math::var;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::algebra_solver_vari;

  int nParms = 2, nSolutions = 2;
  Eigen::VectorXd x(nSolutions);
  x << 1, 1;  // initial guess
  Matrix<var, Dynamic, 1> parms(nParms);
  parms << 36, 6;
  std::vector<double> dummy_dat(0);
  std::vector<int> dummy_dat_int(0);

  algebra_solver_vari<algebraEq_functor<double, double>, var>
   baseVari(algebraEq_functor<double, double>(),
         x, parms, dummy_dat, dummy_dat_int);

  Matrix<double, Dynamic, Dynamic> J(nSolutions, nParms);
  J << 1, 0, 0, 1;

  for (int i = 0; i < nSolutions; i++) {
    stan::math::set_zero_all_adjoints();
    baseVari.impl_->theta_(i).vi_->adj_ = 1;
    baseVari.chain();

    for (int j = 0; j < nParms; j++) {
      EXPECT_EQ(J(i, j), baseVari.impl_->parms_(j).vi_->adj_);
    }
  }
}

template <typename T1, typename T2>
inline Eigen::Matrix<typename boost::math::tools::promote_args<T1, T2>::type,
                     Eigen::Dynamic, 1>
algebraEq2(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
           const Eigen::Matrix<T2, Eigen::Dynamic, 1>& parms,
           const std::vector<double>& dat,
           const std::vector<int>& dat_int) {
  typedef typename boost::math::tools::promote_args<T1, T2>::type scalar;
  Eigen::Matrix<scalar, Eigen::Dynamic, 1> y(3);
  y(0) = x(0) * x(1) * x(2) - parms(0);
  y(1) = x(0) - parms(1);
  y(2) = x(1) - parms(2);
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

TEST(MathMatrix, algebra_solver_chain_2) {
  using stan::math::var;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::algebra_solver_vari;

  int nParms = 3, nSolutions = 3;
  Eigen::VectorXd x(nSolutions);
  x << 1, 1, 1;  // initial guess
  Matrix<var, Dynamic, 1> parms(nParms);
  parms << 5, 4, 2;
  std::vector<double> dummy_dat(0);
  std::vector<int> dummy_dat_int(0);

  algebra_solver_vari<algebraEq2_functor<double, double>, var>
   baseVari(algebraEq2_functor<double, double>(),
         x, parms, dummy_dat, dummy_dat_int);

  double J31 = 1 / (parms(1).val() * parms(2).val());
  double J32 = - parms(0).val() / (parms(1).val() * parms(1).val() * parms(2).val());
  double J33 = - parms(0).val() / (parms(1).val() * parms(2).val() * parms(2).val());

  Matrix<double, Dynamic, Dynamic> J(nSolutions, nParms);
  J << 0, 1, 0, 0, 0, 1, J31, J32, J33;

  for (int i = 0; i < nSolutions; i++) {
    stan::math::set_zero_all_adjoints();
    baseVari.impl_->theta_(i).vi_->adj_ = 1;
    baseVari.chain();

    for (int j = 0; j < nParms; j++) {
      EXPECT_EQ(J(i, j), baseVari.impl_->parms_(j).vi_->adj_);
    }
  }
  } */

template <typename T1, typename T2>
inline Eigen::Matrix<typename boost::math::tools::promote_args<T1, T2>::type,
                     Eigen::Dynamic, 1>
algebraEq3(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
           const Eigen::Matrix<T2, Eigen::Dynamic, 1>& parms,
           const std::vector<double>& dat,
           const std::vector<int>& dat_int) {
  typedef typename boost::math::tools::promote_args<T1, T2>::type scalar;
  Eigen::Matrix<scalar, Eigen::Dynamic, 1> y(2);
  y(0) = x(0) - parms(0) * parms(1); //  + dat[0];
  // y(1) = x(1) - parms(2); //  + dat_int[0] * dat[1];
  y(1) = x(1) - parms(1);
  return y;
}

template <typename T1, typename T2>
struct algebraEq3_functor {
private:
  Eigen::Matrix<T1, Eigen::Dynamic, 1> theta;
  Eigen::Matrix<T2, Eigen::Dynamic, 1> parms;
  std::vector<double> dat;
  std::vector<int> dat_int;
  bool x_is_dv;

public:
  algebraEq3_functor() { };

  algebraEq3_functor(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta_para,
                     const Eigen::Matrix<T2, Eigen::Dynamic, 1>& parms_para,
                     const std::vector<double>& dat_para,
                     const std::vector<int>& dat_int_para,
                     const bool& x_is_dv_para) {
    theta = theta_para;
    parms = parms_para;
    dat = dat_para;
    dat_int = dat_int_para;
    x_is_dv = x_is_dv_para;
  }

  template <typename T>
  inline
  Eigen::Matrix<typename boost::math::tools::promote_args<T, T1, T2>::type,
                Eigen::Dynamic, 1>
  operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1> x) const {
    if (x_is_dv)
      return algebraEq3(x, parms, dat, dat_int);
    else
      return algebraEq3(theta, x, dat, dat_int);
  }
};

TEST(MathMatrix, algebra_solver_chain_3) {
  using stan::math::var;
  using stan::math::vari;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::algebra_solver;
  using stan::math::algebra_solver_vari;
  using stan::math::ChainableStack;
  using stan::math::empty_nested;
  using stan::math::nested_size;

  int nParms = 3, nSolutions = 2;
  Eigen::VectorXd x(nSolutions);
  x << 1, 1;  // initial guess
  Matrix<var, Dynamic, 1> parms(nParms);
  parms << 5, 4, 2;
  std::vector<double> dummy_dat(0);
  std::vector<int> dummy_dat_int(0);

  Matrix<double, Dynamic, Dynamic> J(nSolutions, nParms);
  J << 4, 5, 0, 0, 0, 1;

  Matrix<var, Dynamic, 1> theta
    = algebra_solver(algebraEq3_functor<double, double>(),
                     x, parms, dummy_dat, dummy_dat_int);

  for (int i = 0; i < nSolutions; i++) {
    stan::math::set_zero_all_adjoints();

    // initialize adjoint of final expression
    theta(i).vi_->init_dependent();

    typedef std::vector<vari*>::reverse_iterator it_t;
    it_t begin = ChainableStack::var_stack_.rbegin();
    it_t end = empty_nested()
      ? ChainableStack::var_stack_.rend() : begin + nested_size();  

    std::cout << "Marker A at pass: " << i << std::endl;
    for (it_t it = begin; it < end; ++it) {
      std::cout << "parms_adj before: " 
                << parms(0).vi_->adj_ << " "
                << parms(1).vi_->adj_ << " "
                << parms(2).vi_->adj_ << std::endl;
      (*it)->chain();
      std::cout << "parms_adj after_: " 
                << parms(0).vi_->adj_ << " "
                << parms(1).vi_->adj_ << " "
                << parms(2).vi_->adj_ << std::endl;
    }
    std::cout << "Marker B at pass: " << i << std::endl;

    for (int j = 0; j < nParms; j++) {
      EXPECT_EQ(J(i, j), parms(j).vi_->adj_);
    }

  }
}
