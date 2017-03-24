#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/functor/algebra_solver.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
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
  z(0) = x(0) - y(0);
  z(1) = x(1) - y(1);
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


TEST(MathMatrix, algebra_solver_eq1) {
  using stan::math::algebra_solver;
  using stan::math::sum;
  using stan::math::var;
  using std::cout;
  using std::endl;

  int n_x = 2, n_y = 2;
  for (int k = 0; k < n_x; k++) {
    Eigen::VectorXd x(n_x);
    x << 1, 1;  // initial guess

    Eigen::Matrix<var, Eigen::Dynamic, 1> y(n_x);
    y << 36, 6;

    std::vector<double> dummy_dat(0);
    std::vector<int> dummy_dat_int(0);

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta;
  
    theta = algebra_solver(algebraEq_functor<double, double>(),
                           x, y, dummy_dat, dummy_dat_int);

    EXPECT_EQ(36, theta(0));
    EXPECT_EQ(6, theta(1));
    
    Eigen::MatrixXd J(n_x, n_y);
    J << 1, 0, 0, 1;

    AVEC y_vec = createAVEC(y(0), y(1));
    VEC g;
    theta(k).grad(y_vec, g);

    for (int i = 0; i < n_y; i++)
      EXPECT_EQ(J(k, i), g[i]);
  }
}

///////////////////////////////////////////////////////////////////////////////

template <typename T1, typename T2>
inline Eigen::Matrix<typename boost::math::tools::promote_args<T1, T2>::type,
                     Eigen::Dynamic, 1>
algebraEq2(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
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
struct algebraEq2_functor {
private:
  Eigen::Matrix<T1, Eigen::Dynamic, 1> x_;
  Eigen::Matrix<T2, Eigen::Dynamic, 1> y_;
  std::vector<double> dat_;
  std::vector<int> dat_int_;
  bool x_is_dv_;

public:
  algebraEq2_functor() { };

  algebraEq2_functor(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& x,
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
      return algebraEq2(x, y_, dat_, dat_int_);
    else
      return algebraEq2(x_, x, dat_, dat_int_);
  }
};

TEST(MathMatrix, algebra_solver_eq2) {
  using stan::math::algebra_solver;
  using stan::math::sum;
  using stan::math::var;
  using std::cout;
  using std::endl;

  int n_x = 2, n_y = 3;
  for (int k = 0; k < n_x; k++) {
    Eigen::VectorXd x(n_x);
    x << 1, 1;  // initial guess

    Eigen::Matrix<var, Eigen::Dynamic, 1> y(n_y);
    y << 5, 4, 2;

    std::vector<double> dummy_dat(0);
    std::vector<int> dummy_dat_int(0);

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta;
  
    theta = algebra_solver(algebraEq2_functor<double, double>(),
                           x, y, dummy_dat, dummy_dat_int);

    EXPECT_EQ(20, theta(0));
    EXPECT_EQ(2, theta(1));
    
    Eigen::MatrixXd J(n_x, n_y);
    J << 4, 5, 0, 0, 0, 1;

    AVEC y_vec = createAVEC(y(0), y(1), y(2));
    VEC g;
    theta(k).grad(y_vec, g);

    for (int i = 0; i < n_y; i++)
      EXPECT_EQ(J(k, i), g[i]);
  }
}
