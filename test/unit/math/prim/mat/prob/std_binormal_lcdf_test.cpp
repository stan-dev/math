#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>
#include <iostream>

using Eigen::Dynamic;
using Eigen::Matrix;
using std::vector;

template <typename T>
void to_function_input(int N_y,
                       const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec,
                       vector<Eigen::Matrix<T, Eigen::Dynamic, 1>>& y,
                       vector<T>& rho) {
  int cntr = 0;
  for (int i = 0; i < N_y; ++i) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> el(2);
    el(0) = inp_vec(cntr);
    ++cntr;
    el(1) = inp_vec(cntr);
    ++cntr;
    y.push_back(el);
  }
  for (int i = cntr; i < inp_vec.size(); ++i) {
    rho.push_back(inp_vec(i));
  }
}

template <typename T>
void to_function_input(int N_y,
                       const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec,
                       vector<Eigen::Matrix<T, Eigen::Dynamic, 1>>& y, T& rho) {
  int cntr = 0;
  for (int i = 0; i < N_y; ++i) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> el(2);
    el(0) = inp_vec(cntr);
    ++cntr;
    el(1) = inp_vec(cntr);
    ++cntr;
    y.push_back(el);
  }
  rho = inp_vec(N_y * 2);
}

struct bivar_norm_lcdf {
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    Matrix<T, Dynamic, 1> y(2);
    y(0) = inp_vec(0);
    y(1) = inp_vec(1);
    return stan::math::std_binormal_lcdf(y, inp_vec(2));
  }
};

struct dual_std_vec_bivar_norm_lcdf {
  int N_y_;
  explicit dual_std_vec_bivar_norm_lcdf(int N_y) : N_y_(N_y) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    vector<Matrix<T, Dynamic, 1>> y;
    vector<T> rho;
    to_function_input(N_y_, inp_vec, y, rho);
    return stan::math::std_binormal_lcdf(y, rho);
  }
};

struct one_std_vec_bivar_norm_lcdf {
  int N_y_;
  explicit one_std_vec_bivar_norm_lcdf(int N_y) : N_y_(N_y) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    vector<Matrix<T, Dynamic, 1>> y;
    T rho;
    to_function_input(N_y_, inp_vec, y, rho);
    return stan::math::std_binormal_lcdf(y, rho);
  }
};

struct log_binorm {
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    using std::log;
    return log(
        stan::math::std_binormal_integral(inp_vec(0), inp_vec(1), inp_vec(2)));
  }
};

struct dual_std_vec_log_binorm {
  size_t N_y_;
  explicit dual_std_vec_log_binorm(size_t N_y) : N_y_(N_y) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    using std::log;
    vector<Matrix<T, Dynamic, 1>> y;
    vector<T> rho;
    to_function_input(N_y_, inp_vec, y, rho);
    T accum(0.0);
    if (rho.size() == N_y_) {
      for (size_t i = 0; i < N_y_; ++i)
        accum
            += log(stan::math::std_binormal_integral(y[i][0], y[i][1], rho[i]));
    } else {
      for (size_t i = 0; i < N_y_; ++i)
        accum
            += log(stan::math::std_binormal_integral(y[i][0], y[i][1], rho[0]));
    }
    return accum;
  }
};

TEST(MathFunctions, binormal_lcdf_using) {
  using stan::math::std_binormal_lcdf;
}

TEST(MathFunctions, binormal_integral_throw_RV_1_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  Matrix<double, Dynamic, 1> y(2);
  y << nan, -2.0;
  double rho = 0.3;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_RV_2_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  Matrix<double, Dynamic, 1> y(2);
  y << -2, nan;
  double rho = 0.3;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  Matrix<double, Dynamic, 1> y(2);
  y << -2, 3;
  double rho = nan;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_lt_n_1) {
  Matrix<double, Dynamic, 1> y(2);
  y << -2, 3;
  double rho = -1.3;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_gt_1) {
  Matrix<double, Dynamic, 1> y(2);
  y << -2, 3;
  double rho = 1.3;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_size_y_gt_two) {
  Matrix<double, Dynamic, 1> y(3);
  double rho = 0.3;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::invalid_argument);
}
TEST(MathFunctions, vec_binormal_integral_throw_size_y_gt_two) {
  vector<Matrix<double, Dynamic, 1>> y(2);
  Matrix<double, Dynamic, 1> vec_test(2);
  Matrix<double, Dynamic, 1> vec_test_alt(3);
  y[0] = vec_test;
  y[1] = vec_test_alt;
  double rho = 0.3;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::invalid_argument);
}
TEST(MathFunctions, vec_binormal_integral_rho_y_inconsistent_size) {
  vector<Matrix<double, Dynamic, 1>> y(2);
  vector<double> rho(3);
  rho.push_back(0.4);
  rho.push_back(0.4);
  rho.push_back(0.4);
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::invalid_argument);
}
TEST(MathFunctions, vec_binormal_integral_rho_y_inconsistent_size_eigen) {
  vector<Matrix<double, Dynamic, 1>> y(2);
  Matrix<double, Dynamic, 1> rho(3);
  rho << 0.5, 0.1, 0.2;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::invalid_argument);
}
TEST(MathFunctions, vec_binormal_integral_throw_nan_second_arg) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  vector<Matrix<double, Dynamic, 1>> y(2);
  Matrix<double, Dynamic, 1> vec_test(2);
  Matrix<double, Dynamic, 1> vec_test_alt(2);
  vec_test << 2, nan;
  vec_test_alt << 2, 3;
  y[0] = vec_test;
  y[1] = vec_test_alt;
  double rho = 0.3;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, vec_binormal_integral_throw_nan_second_arg_second_vec) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  vector<Matrix<double, Dynamic, 1>> y(2);
  Matrix<double, Dynamic, 1> vec_test(2);
  Matrix<double, Dynamic, 1> vec_test_alt(2);
  vec_test << 2, nan;
  vec_test_alt << 2, 3;
  y[0] = vec_test_alt;
  y[1] = vec_test;
  double rho = 0.3;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, vec_binormal_integral_throw_nan_first_arg) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  vector<Matrix<double, Dynamic, 1>> y(2);
  Matrix<double, Dynamic, 1> vec_test(2);
  Matrix<double, Dynamic, 1> vec_test_alt(2);
  vec_test << nan, 2;
  vec_test_alt << 2, 3;
  y[0] = vec_test;
  y[1] = vec_test_alt;
  double rho = 0.3;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, vec_binormal_integral_throw_nan_first_arg_second_vec) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  vector<Matrix<double, Dynamic, 1>> y(2);
  Matrix<double, Dynamic, 1> vec_test(2);
  Matrix<double, Dynamic, 1> vec_test_alt(2);
  vec_test << nan, 2;
  vec_test_alt << 2, 3;
  y[0] = vec_test_alt;
  y[1] = vec_test;
  double rho = 0.3;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_no_throw) {
  Matrix<double, Dynamic, 1> y(2);
  y << -2, 3;
  double rho = 0.3;
  EXPECT_NO_THROW(stan::math::std_binormal_lcdf(y, rho));
}
TEST(MathFunctions, binormal_integral_no_throw_vec) {
  vector<Matrix<double, Dynamic, 1>> y(2);
  Matrix<double, Dynamic, 1> y_1(2);
  y_1 << -2, 3;
  y[0] = y_1;
  y[1] = y_1;
  double rho = 0.3;
  EXPECT_NO_THROW(stan::math::std_binormal_lcdf(y, rho));
}
TEST(MathFunctions, binormal_integral_no_throw_y_vec_rho_vec) {
  vector<Matrix<double, Dynamic, 1>> y(2);
  Matrix<double, Dynamic, 1> y_1(2);
  y_1 << -2, 3;
  y[0] = y_1;
  y[1] = y_1;
  Matrix<double, 1, Dynamic> rho(2);
  rho << 0.3, 0.4;
  EXPECT_NO_THROW(stan::math::std_binormal_lcdf(y, rho));
}
TEST(MathFunctions, binormal_integral_no_throw_row_vec) {
  Matrix<double, 1, Dynamic> y(2);
  y << -2, 3;
  double rho = 0.3;
  EXPECT_NO_THROW(stan::math::std_binormal_lcdf(y, rho));
}
TEST(MathFunctions, vec_binormal_integral_no_throw_row_vec) {
  vector<Matrix<double, 1, Dynamic>> y(2);
  Matrix<double, 1, Dynamic> y_1(2);
  y_1 << -2, 3;
  y[0] = y_1;
  y[1] = y_1;
  double rho = 0.3;
  EXPECT_NO_THROW(stan::math::std_binormal_lcdf(y, rho));
}
TEST(MathFunctions, binormal_integral_throw_y_row_vec_rho_vec) {
  Matrix<double, 1, Dynamic> y(2);
  y << -2, 3;
  Matrix<double, 1, Dynamic> rho(4);
  rho << 0.3, 0.4, 0.5, 0.6;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::invalid_argument);
}
TEST(MathFunctions, binormal_integral_throw_vec_y_row_vec_rho_vec) {
  vector<Matrix<double, 1, Dynamic>> y_vec(1);
  Matrix<double, 1, Dynamic> y(2);
  y << -2, 3;
  y_vec[0] = y;
  Matrix<double, 1, Dynamic> rho(4);
  rho << 0.3, 0.4, 0.5, 0.6;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y_vec, rho),
               std::invalid_argument);
}
TEST(MathFunctions, binormal_integral_no_throw_vec_y_row_vec_rho_vec_one_el) {
  vector<Matrix<double, 1, Dynamic>> y_vec(1);
  Matrix<double, 1, Dynamic> y(2);
  y << -2, 3;
  y_vec[0] = y;
  Matrix<double, 1, Dynamic> rho(1);
  rho << 0.3;
  EXPECT_NO_THROW(stan::math::std_binormal_lcdf(y_vec, rho));
}
TEST(MathFunctions,
     binormal_integral_throw_two_el_vec_y_row_vec_rho_vec_one_el) {
  vector<Matrix<double, Dynamic, 1>> y_vec(2);
  Matrix<double, Dynamic, 1> y(2);
  y << -2, 3;
  y_vec[0] = y;
  y_vec[1] = y;
  Matrix<double, 1, Dynamic> rho(1);
  rho << 0.3;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y_vec, rho),
               std::invalid_argument);
}
TEST(MathFunctions, binormal_integral_val_boundaries_test) {
  // Independent normal RVs
  using std::log;
  double rho = 0;
  double a = -0.4;
  double b = 2.7;
  Matrix<double, Dynamic, 1> y(2);
  y(0) = a;
  y(1) = b;
  EXPECT_FLOAT_EQ(log(stan::math::Phi(a)) + log(stan::math::Phi(b)),
                  stan::math::std_binormal_lcdf(y, rho));

  // Perfectly correlated RVs
  rho = 1;
  a = -3.4;
  b = 3.7;
  y(0) = a;
  y(1) = b;
  EXPECT_FLOAT_EQ(log(stan::math::Phi(a)),
                  stan::math::std_binormal_lcdf(y, rho));

  // Perfectly anticorrelated RVs
  rho = -1;
  a = 2.4;
  b = 1.7;
  y(0) = a;
  y(1) = b;
  EXPECT_FLOAT_EQ(log(stan::math::Phi(a) + stan::math::Phi(b) - 1),
                  stan::math::std_binormal_lcdf(y, rho));

  // Perfectly anticorrelated RVs
  rho = -1;
  a = -2.4;
  b = 1.7;
  y(0) = a;
  y(1) = b;
  EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(),
                  stan::math::std_binormal_lcdf(y, rho));

  // a = rho * b
  rho = -0.7;
  b = 1.7;
  a = rho * b;
  y(0) = a;
  y(1) = b;
  EXPECT_FLOAT_EQ(
      log(0.5 / stan::math::pi() * std::exp(-0.5 * b * b) * std::asin(rho)
          + stan::math::Phi(a) * stan::math::Phi(b)),
      stan::math::std_binormal_lcdf(y, rho));

  // b = rho * a
  rho = -0.7;
  a = 1.7;
  b = rho * a;
  y(0) = a;
  y(1) = b;
  EXPECT_FLOAT_EQ(
      log(0.5 / stan::math::pi() * std::exp(-0.5 * a * a) * std::asin(rho)
          + stan::math::Phi(a) * stan::math::Phi(b)),
      stan::math::std_binormal_lcdf(y, rho));
  rho = 0.7;
  a = std::numeric_limits<double>::infinity();
  b = std::numeric_limits<double>::infinity();
  y(0) = a;
  y(1) = b;
  EXPECT_FLOAT_EQ(0, stan::math::std_binormal_lcdf(y, rho));

  rho = 0.7;
  a = -std::numeric_limits<double>::infinity();
  b = -std::numeric_limits<double>::infinity();
  y(0) = a;
  y(1) = b;
  EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(),
                  stan::math::std_binormal_lcdf(y, rho));

  rho = -0.7;
  a = -std::numeric_limits<double>::infinity();
  b = -std::numeric_limits<double>::infinity();
  y(0) = a;
  y(1) = b;
  EXPECT_FLOAT_EQ(a, stan::math::std_binormal_lcdf(y, rho));

  rho = -0.7;
  a = std::numeric_limits<double>::infinity();
  b = std::numeric_limits<double>::infinity();
  y(0) = a;
  y(1) = b;
  EXPECT_FLOAT_EQ(0, stan::math::std_binormal_lcdf(y, rho));

  rho = -0.7;
  a = 1.5;
  b = std::numeric_limits<double>::infinity();
  y(0) = a;
  y(1) = b;
  EXPECT_FLOAT_EQ(log(stan::math::Phi(1.5)),
                  stan::math::std_binormal_lcdf(y, rho));

  rho = 0.7;
  a = 1.5;
  b = std::numeric_limits<double>::infinity();
  y(0) = a;
  y(1) = b;
  EXPECT_FLOAT_EQ(log(stan::math::Phi(1.5)),
                  stan::math::std_binormal_lcdf(y, rho));

  rho = 0.7;
  b = 2.5;
  a = std::numeric_limits<double>::infinity();
  y(0) = a;
  y(1) = b;
  EXPECT_FLOAT_EQ(log(stan::math::Phi(2.5)),
                  stan::math::std_binormal_lcdf(y, rho));

  rho = -0.7;
  b = 0.5;
  a = std::numeric_limits<double>::infinity();
  y(0) = a;
  y(1) = b;
  EXPECT_FLOAT_EQ(log(stan::math::Phi(0.5)),
                  stan::math::std_binormal_lcdf(y, rho));
}
TEST(MathFunctions, binormal_integral_val_test) {
  // Hard-coded values calculated in R using pmvnorm(lower = -Inf, upper = c(a,
  // b), corr = matrix(c(1, rho, rho, 1), 2, 2), algorithm = TVPACK(1e-16))
  // Independent normal RVs
  vector<Matrix<double, Dynamic, 1>> vals;
  Matrix<double, Dynamic, 1> inp_vec(3);
  bivar_norm_lcdf dist_fun;
  log_binorm tru_fun;
  // 000
  inp_vec << 0.4, 2.7, 0.3;
  vals.push_back(inp_vec);
  // 100
  inp_vec << -0.4, 2.7, 0.3;
  vals.push_back(inp_vec);
  // 010
  inp_vec << 0.4, -2.7, 0.3;
  vals.push_back(inp_vec);
  // 110
  inp_vec << -0.4, -2.7, 0.3;
  vals.push_back(inp_vec);
  // 001
  inp_vec << 0.4, 2.7, -0.3;
  vals.push_back(inp_vec);
  // 101
  inp_vec << -0.4, 2.7, -0.3;
  vals.push_back(inp_vec);
  // 011
  inp_vec << 0.4, -2.7, -0.3;
  vals.push_back(inp_vec);
  // 111
  inp_vec << -0.4, -2.7, -0.3;
  vals.push_back(inp_vec);

  inp_vec << -0.4, 2.7, 0.3;
  vals.push_back(inp_vec);
  inp_vec << -0.4, 2.7, 0.99;
  vals.push_back(inp_vec);
  inp_vec << 2.5, 2.7, 0.99;
  vals.push_back(inp_vec);
  inp_vec << 3.5, 3.7, 0.99;
  vals.push_back(inp_vec);
  inp_vec << -4.5, 4.7, -0.99;
  vals.push_back(inp_vec);
  inp_vec << -4.5, 10, -0.99;
  vals.push_back(inp_vec);
  for (size_t i = 0; i < vals.size(); ++i)
    EXPECT_FLOAT_EQ(dist_fun(vals[i]), tru_fun(vals[i]));
}
TEST(MathFunctions, vec_binormal_integral_val_test) {
  // Hard-coded values calculated in R using pmvnorm(lower = -Inf, upper = c(a,
  // b), corr = matrix(c(1, rho, rho, 1), 2, 2), algorithm = TVPACK(1e-16))
  // Independent normal RVs
  int N_y = 3;
  dual_std_vec_bivar_norm_lcdf dist_fun(N_y);
  one_std_vec_bivar_norm_lcdf dist_fun2(N_y);
  dual_std_vec_log_binorm tru_fun(N_y);
  Matrix<double, Dynamic, 1> inp_vec(N_y * 3);
  inp_vec << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7, 0.3, 0.4, 0.5;

  EXPECT_FLOAT_EQ(dist_fun(inp_vec), tru_fun(inp_vec));

  Matrix<double, Dynamic, 1> inp_vec2(N_y * 2 + 1);
  inp_vec2 << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7, 0.3;
  EXPECT_FLOAT_EQ(dist_fun2(inp_vec2), tru_fun(inp_vec2));
}
