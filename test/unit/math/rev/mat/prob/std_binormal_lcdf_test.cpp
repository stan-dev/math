#include <stan/math/rev/mat.hpp>
#include <stan/math/prim/mat/prob/std_binormal_lcdf.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <limits>
#include <vector>

using Eigen::Dynamic;
using Eigen::Matrix;
using std::vector;
using stan::math::var;

template <typename T>
void to_function_input(int N_y, const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec, 
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
  for (int i = cntr; i < inp_vec.size(); ++i)
    rho.push_back(T(inp_vec(i)));
}

template <typename T>
void to_function_input(int N_y, const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec, 
                       vector<Eigen::Matrix<T, Eigen::Dynamic, 1>>& y,
                       T& rho) {
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
  dual_std_vec_bivar_norm_lcdf(int N_y) : N_y_(N_y) {}
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
  one_std_vec_bivar_norm_lcdf(int N_y) : N_y_(N_y) {}
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
    return log(stan::math::std_binormal_integral(inp_vec(0), inp_vec(1), inp_vec(2)));
  }
};

struct dual_std_vec_log_binorm {
  int N_y_;
  dual_std_vec_log_binorm(int N_y) : N_y_(N_y) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    using std::log;
    vector<Matrix<T, Dynamic, 1>> y;
    vector<T> rho;
    to_function_input(N_y_, inp_vec, y, rho);
    T accum(0.0);
    if (rho.size() == N_y_)
      for (int i = 0; i < N_y_; ++i) 
        accum += log(stan::math::std_binormal_integral(y[i][0], y[i][1], rho[i]));
    else 
      for (int i = 0; i < N_y_; ++i) 
        accum += log(stan::math::std_binormal_integral(y[i][0], y[i][1], rho[0]));
    return accum;
  }
};

TEST(MathFunctions, binormal_lcdf_using) { using stan::math::std_binormal_lcdf; }

TEST(MathFunctions, binormal_integral_throw_RV_1_nan) {
  var nan = std::numeric_limits<var>::quiet_NaN();
  Matrix<var, Dynamic, 1> y(2);
  y << nan, -2.0;
  var rho = 0.3;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho),
               std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_RV_2_nan) {
  var nan = std::numeric_limits<var>::quiet_NaN();
  Matrix<var, Dynamic, 1> y(2);
  y << -2.0, nan;
  var rho = 0.3;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho),
               std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_nan) {
  var nan = std::numeric_limits<var>::quiet_NaN();
  Matrix<var, Dynamic, 1> y(2);
  y << -2.0, 0.3;
  var rho = nan;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho),
               std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_lt_neg_one) {
  Matrix<var, Dynamic, 1> y(2);
  y << -2.0, 0.3;
  var rho = -1.3;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho),
               std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_gt_neg_one) {
  Matrix<var, Dynamic, 1> y(2);
  y << -2.0, 0.3;
  var rho = 1.3;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho),
               std::domain_error);
}
TEST(MathFunctions, binormal_integral_val_boundaries_test) {
  // Independent normal RVs
  using std::log;
  using std::exp;
  using std::asin;
  var rho = 0;
  var a = -0.4;
  var b = 2.7;
  Matrix<var, Dynamic, 1> y(2);
  y(0) = a;
  y(1) = b;
  var fun_comp = log(stan::math::Phi(a)) + log(stan::math::Phi(b));
  var fun_run = stan::math::std_binormal_lcdf(y, rho);
  EXPECT_FLOAT_EQ(fun_comp.val(), fun_run.val());

  // Perfectly correlated RVs 
  rho = 1;
  a = -3.4;
  b = 3.7;
  y(0) = a;
  y(1) = b;
  fun_comp = log(stan::math::Phi(a));
  fun_run = stan::math::std_binormal_lcdf(y, rho);
  EXPECT_FLOAT_EQ(fun_comp.val(), fun_run.val());
  

  // Perfectly anticorrelated RVs 
  rho = -1;
  a = 2.4;
  b = 1.7;
  y(0) = a;
  y(1) = b;
  fun_comp = log(stan::math::Phi(a) + stan::math::Phi(b) - 1);
  fun_run = stan::math::std_binormal_lcdf(y, rho);
  
  // Perfectly anticorrelated RVs 
  rho = -1;
  a = -2.4;
  b = 1.7;
  y(0) = a;
  y(1) = b;
  fun_run = stan::math::std_binormal_lcdf(y, rho);
  EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(), fun_run.val());
  
  // a = rho * b 
  rho = -0.7;
  b = 1.7;
  a = rho * b;
  y(0) = a;
  y(1) = b;
  fun_comp = log(0.5 / stan::math::pi()
                  * exp(-0.5 * b * b)
                  * asin(rho)
                  + stan::math::Phi(a) * stan::math::Phi(b));
  fun_run = stan::math::std_binormal_lcdf(y, rho);
  EXPECT_FLOAT_EQ(fun_comp.val(), fun_run.val());
  
  // b = rho * a 
  rho = -0.7;
  a = 1.7;
  b = rho * a;
  y(0) = a;
  y(1) = b;
  fun_comp = log(0.5 / stan::math::pi()
                  * exp(-0.5 * a * a)
                  * asin(rho)
                  + stan::math::Phi(a) * stan::math::Phi(b));
  fun_run = stan::math::std_binormal_lcdf(y, rho);
  EXPECT_FLOAT_EQ(fun_comp.val(), fun_run.val());

  rho = 0.7;
  a = std::numeric_limits<double>::infinity();
  b = std::numeric_limits<double>::infinity();
  y(0) = a;
  y(1) = b;
  fun_run = stan::math::std_binormal_lcdf(y, rho);
  EXPECT_FLOAT_EQ(0,fun_run.val()); 

  rho = 0.7;
  a = -std::numeric_limits<double>::infinity();
  b = -std::numeric_limits<double>::infinity();
  y(0) = a;
  y(1) = b;
  fun_run = stan::math::std_binormal_lcdf(y, rho);
  EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(),fun_run.val()); 

  rho = -0.7;
  a = -std::numeric_limits<double>::infinity();
  b = -std::numeric_limits<double>::infinity();
  y(0) = a;
  y(1) = b;
  fun_run = stan::math::std_binormal_lcdf(y, rho);
  EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(),fun_run.val()); 

  rho = -0.7;
  a = std::numeric_limits<double>::infinity();
  b = std::numeric_limits<double>::infinity();
  y(0) = a;
  y(1) = b;
  fun_run = stan::math::std_binormal_lcdf(y, rho);
  EXPECT_FLOAT_EQ(0,fun_run.val()); 

  rho = -0.7;
  a = 1.5;
  b = std::numeric_limits<double>::infinity();
  y(0) = a;
  y(1) = b;
  fun_run = stan::math::std_binormal_lcdf(y, rho);
  fun_comp = log(stan::math::Phi(1.5));
  EXPECT_FLOAT_EQ(fun_comp.val(), fun_run.val());

  rho = 0.7;
  a = 1.5;
  b = std::numeric_limits<double>::infinity();
  y(0) = a;
  y(1) = b;
  fun_comp = log(stan::math::Phi(1.5));
  fun_run = stan::math::std_binormal_lcdf(y, rho);
  EXPECT_FLOAT_EQ(fun_comp.val(), fun_run.val());

  rho = 0.7;
  b = 2.5;
  a = std::numeric_limits<double>::infinity();
  y(0) = a;
  y(1) = b;
  fun_comp = log(stan::math::Phi(2.5));
  fun_run = stan::math::std_binormal_lcdf(y, rho);
  EXPECT_FLOAT_EQ(fun_comp.val(), fun_run.val());

  rho = -0.7;
  b = 0.5;
  a = std::numeric_limits<double>::infinity();
  y(0) = a;
  y(1) = b;
  fun_comp = log(stan::math::Phi(0.5));
  fun_run = stan::math::std_binormal_lcdf(y, rho);
  EXPECT_FLOAT_EQ(fun_comp.val(), fun_run.val());
}
TEST(MathFunctions, binormal_integral_val_test) {
  // Hard-coded values calculated in R using pmvnorm(lower = -Inf, upper = c(a,b), 
  // corr = matrix(c(1,rho,rho,1),2,2), algorithm = TVPACK(1e-16))
  // Independent normal RVs
  vector<Matrix<var, Dynamic, 1>> vals(14, Matrix<var, Dynamic, 1>(3));
  bivar_norm_lcdf dist_fun;
  log_binorm tru_fun;
  // 000
  vals[0] << 0.4, 2.7, 0.3; 
  // 100
  vals[1] << -0.4, 2.7, 0.3; 
  // 010
  vals[2] << 0.4, -2.7, 0.3; 
  // 110
  vals[3] << -0.4, -2.7, 0.3; 
  // 001
  vals[4] << 0.4, 2.7, -0.3; 
  // 101
  vals[5] << -0.4, 2.7, -0.3; 
  // 011
  vals[6] << 0.4, -2.7, -0.3; 
  // 111
  vals[7] << -0.4, -2.7, -0.3; 

  vals[8] << -0.4, 2.7, 0.3; 
  vals[9] << -0.4, 2.7, 0.99; 
  vals[10] << 2.5, 2.7, 0.99; 
  vals[11] << 3.5, 3.7, 0.99; 
  vals[12] << -4.5, 4.7, -0.99; 
  vals[13] << -4.5, 10, -0.99; 
  for (size_t i = 0; i < vals.size(); ++i)
    EXPECT_FLOAT_EQ(dist_fun(vals[i]).val(), tru_fun(vals[i]).val());
}
TEST(MathFunctions, vec_binormal_integral_val_test) {
  // Hard-coded values calculated in R using pmvnorm(lower = -Inf, upper = c(a,b), 
  // corr = matrix(c(1,rho,rho,1),2,2), algorithm = TVPACK(1e-16))
  // Independent normal RVs
  int N_y = 3;
  dual_std_vec_bivar_norm_lcdf dist_fun(N_y);
  one_std_vec_bivar_norm_lcdf dist_fun2(N_y);
  dual_std_vec_log_binorm tru_fun(N_y);
  Matrix<var, Dynamic, 1> inp_vec(N_y * 3);
  inp_vec << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7, 0.3, 0.4, 0.5;

  EXPECT_FLOAT_EQ(dist_fun(inp_vec).val(), tru_fun(inp_vec).val());

  Matrix<var, Dynamic, 1> inp_vec2(N_y * 2 + 1);
  inp_vec2 << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7, 0.3;
  EXPECT_FLOAT_EQ(dist_fun2(inp_vec2).val(), tru_fun(inp_vec2).val());
}
TEST(MathFunctions, vec_binormal_integral_grad_test) {
  // Hard-coded values calculated in R using pmvnorm(lower = -Inf, upper = c(a,b), 
  // corr = matrix(c(1,rho,rho,1),2,2), algorithm = TVPACK(1e-16))
  // Independent normal RVs
  int N_y = 3;
  dual_std_vec_bivar_norm_lcdf dist_fun(N_y);
  one_std_vec_bivar_norm_lcdf dist_fun2(N_y);
  dual_std_vec_log_binorm tru_fun(N_y);
  Matrix<double, Dynamic, 1> inp_vec(N_y * 3);
  inp_vec << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7, 0.3, 0.4, 0.5;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun, inp_vec, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec, fx_tru, grad_fx_tru);

  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_grad2_test) {
  // Hard-coded values calculated in R using pmvnorm(lower = -Inf, upper = c(a,b), 
  // corr = matrix(c(1,rho,rho,1),2,2), algorithm = TVPACK(1e-16))
  // Independent normal RVs
  int N_y = 3;
  one_std_vec_bivar_norm_lcdf dist_fun2(N_y);
  dual_std_vec_log_binorm tru_fun(N_y);
  Matrix<double, Dynamic, 1> inp_vec2(N_y * 2 + 1);
  inp_vec2 << 0.4, -2.7,0.4, -2.7,0.4, -2.7, 0.3;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun2, inp_vec2, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec2, fx_tru, grad_fx_tru);

  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i));
  }
}
