#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <iostream>
#include <sstream>
#include <vector>

struct f1 {
  template <typename T1, typename T2>
  inline
  typename stan::return_type<T1, T2>::type
  operator()(const T1& x, const T2& y, std::ostream& msgs) const {
    return exp(x) + y;
  }
};

struct f2 {
  template <typename T1, typename T2>
  inline
  typename stan::return_type<T1, T2>::type
  operator()(const T1& x, const std::vector<T2>& y, std::ostream& msgs)
  const {
    return exp(x) + pow(y[0], 2) + pow(y[1], 3);
  }
};

struct f3 {
  template <typename T1, typename T2>
  inline
  typename stan::return_type<T1, T2>::type
  operator()(const T1& x, const Eigen::Matrix<T2, Eigen::Dynamic, 1>& y,
             std::ostream& msgs) const {
    return exp(x) + pow(y(0), 2) + pow(y(1), 4) + 3*y(2);
  }
};

struct f4 {
  template <typename T1, typename T2>
  inline
  typename stan::return_type<T1, T2>::type
  operator()(const T1& x, const Eigen::Matrix<T2, 1, Eigen::Dynamic>& y,
             std::ostream& msgs) const {
    return exp(x) + pow(y(0), 2) + pow(y(1), 5) + 3*y(2);
  }
};

std::ostringstream msgs;

TEST(StanMath_integrate_1d_tsc, test1) {
  using stan::math::integrate_1d_tsc;

  f1 if1;

  EXPECT_FLOAT_EQ(integrate_1d_tsc(if1, .2, .7, .5, msgs), 1.04235);

  f2 if2;

  EXPECT_FLOAT_EQ(integrate_1d_tsc(if2, -.2, .7,
                                     std::vector<double>(2, .4), msgs),
                  1.396622);
  f3 if3;
  Eigen::VectorXd a(3);
  a(0) = 4.0;
  a(1) = 6.0;
  a(2) = 5.1;

  EXPECT_FLOAT_EQ(integrate_1d_tsc(if3, -.2, 2.9, a, msgs),
                  4131.985);

  f4 if4;
  Eigen::RowVectorXd b(3);
  b(0) = 4.0;
  b(1) = 6.0;
  b(2) = 5.1;

  EXPECT_FLOAT_EQ(integrate_1d_tsc(if4, -.2, 2.9, b, msgs),
                  24219.99);
}
