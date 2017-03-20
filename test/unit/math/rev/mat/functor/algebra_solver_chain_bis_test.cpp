#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/functor/algebra_solver.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <stan/math/rev/core.hpp>
#include <gtest/gtest.h>
#include <iostream>

template <typename T1, typename T2>
inline Eigen::Matrix<typename boost::math::tools::promote_args<T1, T2>::type,
                     Eigen::Dynamic, 1>
algebraEq(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
          const Eigen::Matrix<T2, Eigen::Dynamic, 1>& y,
          const std::vector<double>& dat,
          const std::vector<int>& dat_int) {
  typedef typename boost::math::tools::promote_args<T1, T2>::type scalar;
  Eigen::Matrix<scalar, Eigen::Dynamic, 1> z(2);
  assert(z.rows() == x.rows());
  z(0) = x(0) - y(0) * y(1);
  z(1) = x(1) - y(2);
  return z;
}

template <typename T1, typename T2>
struct algebraEq_functor {
private:
  Eigen::Matrix<T1, Eigen::Dynamic, 1> x_;
  Eigen::Matrix<T2, Eigen::Dynamic, 1> y_;
  std::vector<double> dat_;
  std::vector<int> dat_int_;
  bool x_is_dv_;

public:
  algebraEq_functor() { };

  algebraEq_functor(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
                    const Eigen::Matrix<T2, Eigen::Dynamic, 1>& y,
                    const std::vector<double>& dat,
                    const std::vector<int>& dat_int,
                    const bool& x_is_dv) {
    x_ = x;
    y_ = y;
    dat_ = dat;
    dat_int_ = dat_int;
    x_is_dv_ = x_is_dv;
  }

  template <typename T>
  inline
  Eigen::Matrix<typename boost::math::tools::promote_args<T, T1, T2>::type,
                Eigen::Dynamic, 1>
  operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1> x) const {
    if (x_is_dv_)
      return algebraEq(x, y_, dat_, dat_int_);
    else
      return algebraEq(x_, x, dat_, dat_int_);
  }
};

TEST(MathMatrix, algebra_solver_chain_bis) {
  using stan::math::algebra_solver;
  using stan::math::algebra_solver_vari;
  using stan::math::ChainableStack;
  using stan::math::empty_nested;
  using stan::math::nested_size;
  using stan::math::var;
  using stan::math::vari;
  using Eigen::VectorXd;
  using Eigen::MatrixXd;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using std::vector;

  int nParms = 3, nSolutions = 2;
  VectorXd x(nSolutions);
  x << 1, 1;
  Matrix<var, Dynamic, 1> y(nParms);
  y << 5, 4, 2;
  vector<double> dummy_dat(0);
  vector<int> dummy_dat_int(0);

  MatrixXd J(nSolutions, nParms);
  J << 4, 5, 0, 0, 0, 1;

  Matrix<var, Dynamic, 1> theta
    = algebra_solver(algebraEq_functor<double, double>(),
                     x, y, dummy_dat, dummy_dat_int);

  EXPECT_EQ(20, theta(0));
  EXPECT_EQ(2, theta(1));

  for (int i = 0; i < nSolutions; i++) {
    stan::math::set_zero_all_adjoints();
    theta(i).vi_->init_dependent();

    typedef vector<vari*>::reverse_iterator it_t;
    it_t begin = ChainableStack::var_stack_.rbegin();

    it_t it = begin;  // CHECK: Is this variable well-defined ??
    (*it)->chain();
    ++it;
    std::cout << "Marker A at pass: " << i << std::endl;
    std::cout << it[0] << " ";
    std::cout << it[1] << " ";
    std::cout << it[2] << std::endl;
    std::cout << "Marker B at pass: " << i << std::endl;
    (*it)->chain();  // BUG: this line of code breaks
    std::cout << "Marker C at pass: " << i << std::endl;
    
    // TEST: the above line of code breaking is equivalent to
    // the FOR loop breaking during the second iteration.

  /*
  it_t end = empty_nested()
    ? ChainableStack::var_stack_.rend() : begin + nested_size();

  for (it_t it = begin; it < end; ++it) {
    std::cout << "parms_adj before: " 
              << y(0).vi_->adj_ << " "
              << y(1).vi_->adj_ << " "
              << y(2).vi_->adj_ << std::endl;
    (*it)->chain();
    std::cout << "parms_adj after_: " 
              << y(0).vi_->adj_ << " "
              << y(1).vi_->adj_ << " "
              << y(2).vi_->adj_ << std::endl;
  } */

    for (int j = 0; j < nParms; j++)
      EXPECT_EQ(J(i, j), y(j).vi_->adj_);
  
  }
}  
  
  
  
