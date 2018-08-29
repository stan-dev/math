#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <limits>
#include <vector>

using Eigen::Dynamic;
using Eigen::Matrix;
using stan::math::var;
using std::vector;

// nector, real
template <typename T>
void to_function_input(int N_y,
                       const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec,
                       Eigen::Matrix<T, Eigen::Dynamic, 1>& y, T& rho) {
  y(0) = inp_vec(0);
  y(1) = inp_vec(1);
  rho = T(inp_vec(2));
}

// Row Vector, real
template <typename T>
void to_function_input(int N_y,
                       const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec,
                       Eigen::Matrix<T, 1, Eigen::Dynamic>& y, T& rho) {
  y(0) = inp_vec(0);
  y(1) = inp_vec(1);
  rho = T(inp_vec(2));
}

// vector<Vector>, real
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

// vector<Row Vector>, real
template <typename T>
void to_function_input(int N_y,
                       const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec,
                       vector<Eigen::Matrix<T, 1, Eigen::Dynamic>>& y, T& rho) {
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

// vector<Vector>, vector
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
  for (int i = cntr; i < inp_vec.size(); ++i)
    rho.push_back(T(inp_vec(i)));
}

// vector<Row Vector>, vector
template <typename T>
void to_function_input(int N_y,
                       const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec,
                       vector<Eigen::Matrix<T, 1, Eigen::Dynamic>>& y,
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

// vector<Vector>, Vector
template <typename T>
void to_function_input(int N_y,
                       const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec,
                       vector<Eigen::Matrix<T, Eigen::Dynamic, 1>>& y,
                       Eigen::Matrix<T, Eigen::Dynamic, 1>& rho) {
  int cntr = 0;
  for (int i = 0; i < N_y; ++i) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> el(2);
    el(0) = inp_vec(cntr);
    ++cntr;
    el(1) = inp_vec(cntr);
    ++cntr;
    y.push_back(el);
  }
  rho.resize(inp_vec.size() - 2 * N_y);
  int j = 0;
  for (int i = cntr; i < inp_vec.size(); ++i)
    rho(j++) = (T(inp_vec(i)));
}

// vector<Row Vector>, Vector
template <typename T>
void to_function_input(int N_y,
                       const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec,
                       vector<Eigen::Matrix<T, 1, Eigen::Dynamic>>& y,
                       Eigen::Matrix<T, Eigen::Dynamic, 1>& rho) {
  int cntr = 0;
  for (int i = 0; i < N_y; ++i) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> el(2);
    el(0) = inp_vec(cntr);
    ++cntr;
    el(1) = inp_vec(cntr);
    ++cntr;
    y.push_back(el);
  }
  rho.resize(inp_vec.size() - 2 * N_y);
  int j = 0;
  for (int i = cntr; i < inp_vec.size(); ++i)
    rho(j++) = (T(inp_vec(i)));
}

// Row vector, nada
template <typename T>
void to_function_input(int N_y,
                       const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec,
                       Eigen::Matrix<T, 1, Eigen::Dynamic>& y) {
  y(0) = inp_vec(0);
  y(1) = inp_vec(1);
}

// vector<Vector>, nada
template <typename T>
void to_function_input(int N_y,
                       const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec,
                       vector<Eigen::Matrix<T, Eigen::Dynamic, 1>>& y) {
  int cntr = 0;
  for (int i = 0; i < N_y; ++i) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> el(2);
    el(0) = inp_vec(cntr);
    ++cntr;
    el(1) = inp_vec(cntr);
    ++cntr;
    y.push_back(el);
  }
}

// vector<Row Vector>, nada
template <typename T>
void to_function_input(int N_y,
                       const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec,
                       vector<Eigen::Matrix<T, 1, Eigen::Dynamic>>& y) {
  int cntr = 0;
  for (int i = 0; i < N_y; ++i) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> el(2);
    el(0) = inp_vec(cntr);
    ++cntr;
    el(1) = inp_vec(cntr);
    ++cntr;
    y.push_back(el);
  }
}

// Vector, real
struct V_R_std_binorm_lcdf {
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    Matrix<T, Dynamic, 1> y(2);
    y(0) = inp_vec(0);
    y(1) = inp_vec(1);
    return stan::math::std_binormal_lcdf(y, inp_vec(2));
  }
};

// Row Vector, real
struct RV_R_std_binorm_lcdf {
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    Matrix<T, 1, Dynamic> y(2);
    y(0) = inp_vec(0);
    y(1) = inp_vec(1);
    return stan::math::std_binormal_lcdf(y, inp_vec(2));
  }
};

// vector<Vector>, real
struct vV_R_std_binorm_lcdf {
  int N_y_;
  explicit vV_R_std_binorm_lcdf(int N_y) : N_y_(N_y) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    vector<Matrix<T, Dynamic, 1>> y;
    T rho;
    to_function_input(N_y_, inp_vec, y, rho);
    return stan::math::std_binormal_lcdf(y, rho);
  }
};

// vector<Row Vector>, real
struct vRV_R_std_binorm_lcdf {
  int N_y_;
  explicit vRV_R_std_binorm_lcdf(int N_y) : N_y_(N_y) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    vector<Matrix<T, 1, Dynamic>> y;
    T rho;
    to_function_input(N_y_, inp_vec, y, rho);
    return stan::math::std_binormal_lcdf(y, rho);
  }
};

// vector<Vector>, vector
struct vV_v_std_binorm_lcdf {
  int N_y_;
  explicit vV_v_std_binorm_lcdf(int N_y) : N_y_(N_y) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    vector<Matrix<T, Dynamic, 1>> y;
    vector<T> rho;
    to_function_input(N_y_, inp_vec, y, rho);
    return stan::math::std_binormal_lcdf(y, rho);
  }
};

// vector<Row Vector>, vector
struct vRV_v_std_binorm_lcdf {
  int N_y_;
  explicit vRV_v_std_binorm_lcdf(int N_y) : N_y_(N_y) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    vector<Matrix<T, 1, Dynamic>> y;
    vector<T> rho;
    to_function_input(N_y_, inp_vec, y, rho);
    return stan::math::std_binormal_lcdf(y, rho);
  }
};

// vector<Vector>, Vector
struct vV_V_std_binorm_lcdf {
  int N_y_;
  explicit vV_V_std_binorm_lcdf(int N_y) : N_y_(N_y) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    vector<Matrix<T, Dynamic, 1>> y;
    Matrix<T, Dynamic, 1> rho;
    to_function_input(N_y_, inp_vec, y, rho);
    return stan::math::std_binormal_lcdf(y, rho);
  }
};

// vector<Row Vector>, Vector
struct vRV_V_std_binorm_lcdf {
  int N_y_;
  explicit vRV_V_std_binorm_lcdf(int N_y) : N_y_(N_y) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    vector<Matrix<T, 1, Dynamic>> y;
    Matrix<T, Dynamic, 1> rho;
    to_function_input(N_y_, inp_vec, y, rho);
    return stan::math::std_binormal_lcdf(y, rho);
  }
};

// Vector, double
struct V_D_std_binorm_lcdf {
  double rho_;
  explicit V_D_std_binorm_lcdf(double rho) : rho_(rho) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    Matrix<T, Dynamic, 1> y(2);
    y(0) = inp_vec(0);
    y(1) = inp_vec(1);
    return stan::math::std_binormal_lcdf(y, rho_);
  }
};

// Row Vector, double
struct RV_D_std_binorm_lcdf {
  double rho_;
  explicit RV_D_std_binorm_lcdf(double rho) : rho_(rho) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    Matrix<T, 1, Dynamic> y(2);
    y(0) = inp_vec(0);
    y(1) = inp_vec(1);
    return stan::math::std_binormal_lcdf(y, rho_);
  }
};

// vector<Vector>, double
struct vV_D_std_binorm_lcdf {
  int N_y_;
  double rho_;
  vV_D_std_binorm_lcdf(int N_y, double rho) : N_y_(N_y), rho_(rho) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    vector<Matrix<T, Dynamic, 1>> y;
    to_function_input(N_y_, inp_vec, y);
    return stan::math::std_binormal_lcdf(y, rho_);
  }
};

// vector<Row Vector>, double
struct vRV_D_std_binorm_lcdf {
  int N_y_;
  double rho_;
  vRV_D_std_binorm_lcdf(int N_y, double rho) : N_y_(N_y), rho_(rho) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    vector<Matrix<T, 1, Dynamic>> y;
    to_function_input(N_y_, inp_vec, y);
    return stan::math::std_binormal_lcdf(y, rho_);
  }
};

// vector<Vector>, double vector
struct vV_vD_std_binorm_lcdf {
  int N_y_;
  vector<double> rho_;
  vV_vD_std_binorm_lcdf(int N_y, vector<double>& rho) : N_y_(N_y), rho_(rho) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    vector<Matrix<T, Dynamic, 1>> y;
    to_function_input(N_y_, inp_vec, y);
    return stan::math::std_binormal_lcdf(y, rho_);
  }
};

// vector<Row Vector>, double vector
struct vRV_vD_std_binorm_lcdf {
  int N_y_;
  vector<double> rho_;
  vRV_vD_std_binorm_lcdf(int N_y, vector<double>& rho) : N_y_(N_y), rho_(rho) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    vector<Matrix<T, 1, Dynamic>> y;
    to_function_input(N_y_, inp_vec, y);
    return stan::math::std_binormal_lcdf(y, rho_);
  }
};

// vector<Vector>, Vector
struct vV_VD_std_binorm_lcdf {
  int N_y_;
  Matrix<double, Dynamic, 1> rho_;
  vV_VD_std_binorm_lcdf(int N_y, Matrix<double, Dynamic, 1>& rho)
      : N_y_(N_y), rho_(rho) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    vector<Matrix<T, Dynamic, 1>> y;
    to_function_input(N_y_, inp_vec, y);
    return stan::math::std_binormal_lcdf(y, rho_);
  }
};

// vector<Row Vector>, Vector
struct vRV_VD_std_binorm_lcdf {
  int N_y_;
  Matrix<double, Dynamic, 1> rho_;
  vRV_VD_std_binorm_lcdf(int N_y, Matrix<double, Dynamic, 1>& rho)
      : N_y_(N_y), rho_(rho) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    vector<Matrix<T, 1, Dynamic>> y;
    to_function_input(N_y_, inp_vec, y);
    return stan::math::std_binormal_lcdf(y, rho_);
  }
};

// double Vector, real
struct VD_R_std_binorm_lcdf {
  Matrix<double, Dynamic, 1> y_;
  explicit VD_R_std_binorm_lcdf(Matrix<double, Dynamic, 1>& y) : y_(y) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    return stan::math::std_binormal_lcdf(y_, inp_vec(0));
  }
};

// double Row Vector, real
struct RVD_R_std_binorm_lcdf {
  Matrix<double, 1, Dynamic> y_;
  explicit RVD_R_std_binorm_lcdf(Matrix<double, 1, Dynamic>& y) : y_(y) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    return stan::math::std_binormal_lcdf(y_, inp_vec(0));
  }
};

// double vector<Vector>, real
struct vVD_R_std_binorm_lcdf {
  vector<Matrix<double, Dynamic, 1>> y_;
  explicit vVD_R_std_binorm_lcdf(vector<Matrix<double, Dynamic, 1>>& y) {
    for (size_t i = 0; i < y.size(); ++i)
      y_.push_back(y[i]);
  }
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    vector<Matrix<T, Dynamic, 1>> y;
    return stan::math::std_binormal_lcdf(y_, inp_vec(0));
  }
};

// double vector<Row Vector>, real
struct vRVD_R_std_binorm_lcdf {
  vector<Matrix<double, 1, Dynamic>> y_;
  explicit vRVD_R_std_binorm_lcdf(vector<Matrix<double, 1, Dynamic>>& y) {
    for (size_t i = 0; i < y.size(); ++i)
      y_.push_back(y[i]);
  }
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    return stan::math::std_binormal_lcdf(y_, inp_vec(0));
  }
};

// double vector<Vector>, vector
struct vVD_v_std_binorm_lcdf {
  vector<Matrix<double, Dynamic, 1>> y_;
  explicit vVD_v_std_binorm_lcdf(vector<Matrix<double, Dynamic, 1>>& y) {
    for (size_t i = 0; i < y.size(); ++i)
      y_.push_back(y[i]);
  }
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    vector<T> rho;
    for (int i = 0; i < inp_vec.size(); ++i)
      rho.push_back(inp_vec(i));
    return stan::math::std_binormal_lcdf(y_, rho);
  }
};

// double vector<Row Vector>, vector
struct vRVD_v_std_binorm_lcdf {
  vector<Matrix<double, 1, Dynamic>> y_;
  explicit vRVD_v_std_binorm_lcdf(vector<Matrix<double, 1, Dynamic>>& y) {
    for (size_t i = 0; i < y.size(); ++i)
      y_.push_back(y[i]);
  }
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    vector<T> rho;
    for (int i = 0; i < inp_vec.size(); ++i)
      rho.push_back(inp_vec(i));
    return stan::math::std_binormal_lcdf(y_, rho);
  }
};

// double vector<Row Vector>, Vector
struct vRVD_V_std_binorm_lcdf {
  vector<Matrix<double, 1, Dynamic>> y_;
  explicit vRVD_V_std_binorm_lcdf(vector<Matrix<double, 1, Dynamic>>& y) {
    for (size_t i = 0; i < y.size(); ++i)
      y_.push_back(y[i]);
  }
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    return stan::math::std_binormal_lcdf(y_, inp_vec);
  }
};

// double vector<Vector>, Vector
struct vVD_V_std_binorm_lcdf {
  vector<Matrix<double, Dynamic, 1>> y_;
  explicit vVD_V_std_binorm_lcdf(vector<Matrix<double, Dynamic, 1>>& y) {
    for (size_t i = 0; i < y.size(); ++i)
      y_.push_back(y[i]);
  }
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    return stan::math::std_binormal_lcdf(y_, inp_vec);
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

struct vV_vD_D_log_std_binorm_integral {
  size_t N_y_;
  vector<double> rho_;
  vV_vD_D_log_std_binorm_integral(size_t N_y, vector<double>& rho)
      : N_y_(N_y), rho_(rho) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    using std::log;
    vector<Matrix<T, Dynamic, 1>> y;
    to_function_input(N_y_, inp_vec, y);
    T accum(0.0);
    if (rho_.size() == N_y_) {
      for (size_t i = 0; i < N_y_; ++i)
        accum += log(
            stan::math::std_binormal_integral(y[i][0], y[i][1], rho_[i]));
    } else {
      for (size_t i = 0; i < N_y_; ++i)
        accum += log(
            stan::math::std_binormal_integral(y[i][0], y[i][1], rho_[0]));
    }
    return accum;
  }
};

struct vV_v_R_log_std_binorm_integral {
  size_t N_y_;
  explicit vV_v_R_log_std_binorm_integral(size_t N_y) : N_y_(N_y) {}
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

template <int R, int C>
struct vVD_v_log_std_binorm_integral {
  vector<Matrix<double, R, C>> y_;
  explicit vVD_v_log_std_binorm_integral(vector<Matrix<double, R, C>>& y) {
    for (size_t i = 0; i < y.size(); ++i)
      y_.push_back(y[i]);
  }
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    using std::log;
    vector<T> rho;
    for (int i = 0; i < inp_vec.size(); ++i)
      rho.push_back(inp_vec(i));
    T accum(0.0);
    if (rho.size() == y_.size()) {
      for (size_t i = 0; i < y_.size(); ++i)
        accum += log(
            stan::math::std_binormal_integral(y_[i][0], y_[i][1], rho[i]));
    } else {
      for (size_t i = 0; i < y_.size(); ++i)
        accum += log(
            stan::math::std_binormal_integral(y_[i][0], y_[i][1], rho[0]));
    }
    return accum;
  }
};

TEST(MathFunctions, binormal_lcdf_using) {
  using stan::math::std_binormal_lcdf;
}

TEST(MathFunctions, binormal_integral_throw_RV_1_nan) {
  var nan = std::numeric_limits<var>::quiet_NaN();
  Matrix<var, Dynamic, 1> y(2);
  y << nan, -2.0;
  var rho = 0.3;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_RV_2_nan) {
  var nan = std::numeric_limits<var>::quiet_NaN();
  Matrix<var, Dynamic, 1> y(2);
  y << -2.0, nan;
  var rho = 0.3;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_nan) {
  var nan = std::numeric_limits<var>::quiet_NaN();
  Matrix<var, Dynamic, 1> y(2);
  y << -2.0, 0.3;
  var rho = nan;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_lt_neg_one) {
  Matrix<var, Dynamic, 1> y(2);
  y << -2.0, 0.3;
  var rho = -1.3;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_gt_neg_one) {
  Matrix<var, Dynamic, 1> y(2);
  y << -2.0, 0.3;
  var rho = 1.3;
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_val_boundaries_test) {
  // Independent normal RVs
  using std::asin;
  using std::exp;
  using std::log;
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
  fun_comp = log(0.5 / stan::math::pi() * exp(-0.5 * b * b) * asin(rho)
                 + stan::math::Phi(a) * stan::math::Phi(b));
  fun_run = stan::math::std_binormal_lcdf(y, rho);
  EXPECT_FLOAT_EQ(fun_comp.val(), fun_run.val());

  // b = rho * a
  rho = -0.7;
  a = 1.7;
  b = rho * a;
  y(0) = a;
  y(1) = b;
  fun_comp = log(0.5 / stan::math::pi() * exp(-0.5 * a * a) * asin(rho)
                 + stan::math::Phi(a) * stan::math::Phi(b));
  fun_run = stan::math::std_binormal_lcdf(y, rho);
  EXPECT_FLOAT_EQ(fun_comp.val(), fun_run.val());

  rho = 0.7;
  a = std::numeric_limits<double>::infinity();
  b = std::numeric_limits<double>::infinity();
  y(0) = a;
  y(1) = b;
  fun_run = stan::math::std_binormal_lcdf(y, rho);
  EXPECT_FLOAT_EQ(0, fun_run.val());

  rho = 0.7;
  a = -std::numeric_limits<double>::infinity();
  b = -std::numeric_limits<double>::infinity();
  y(0) = a;
  y(1) = b;
  fun_run = stan::math::std_binormal_lcdf(y, rho);
  EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(), fun_run.val());

  rho = -0.7;
  a = -std::numeric_limits<double>::infinity();
  b = -std::numeric_limits<double>::infinity();
  y(0) = a;
  y(1) = b;
  fun_run = stan::math::std_binormal_lcdf(y, rho);
  EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(), fun_run.val());

  rho = -0.7;
  a = std::numeric_limits<double>::infinity();
  b = std::numeric_limits<double>::infinity();
  y(0) = a;
  y(1) = b;
  fun_run = stan::math::std_binormal_lcdf(y, rho);
  EXPECT_FLOAT_EQ(0, fun_run.val());

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
  vector<Matrix<var, Dynamic, 1>> vals(14, Matrix<var, Dynamic, 1>(3));
  V_R_std_binorm_lcdf dist_fun;
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
TEST(MathFunctions, vec_binormal_integral_val_test_vV_v) {
  int N_y = 3;
  vV_v_std_binorm_lcdf dist_fun(N_y);
  vV_R_std_binorm_lcdf dist_fun2(N_y);
  vV_v_R_log_std_binorm_integral tru_fun(N_y);
  Matrix<var, Dynamic, 1> inp_vec(N_y * 3);
  inp_vec << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7, 0.3, 0.4, 0.5;

  EXPECT_FLOAT_EQ(dist_fun(inp_vec).val(), tru_fun(inp_vec).val());

  Matrix<var, Dynamic, 1> inp_vec2(N_y * 2 + 1);
  inp_vec2 << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7, 0.3;
  EXPECT_FLOAT_EQ(dist_fun2(inp_vec2).val(), tru_fun(inp_vec2).val());
}
TEST(MathFunctions, vec_binormal_integral_grad_test_V_R) {
  V_R_std_binorm_lcdf dist_fun;
  log_binorm tru_fun;
  Matrix<double, Dynamic, 1> inp_vec(3);
  inp_vec << 0.4, -2.7, 0.8;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun, inp_vec, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun, inp_vec, fx_fd, grad_fx_fd);

  EXPECT_FLOAT_EQ(fx_test, fx_tru);
  EXPECT_FLOAT_EQ(fx_test, fx_fd);
  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i));
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_grad_test_RV_R) {
  RV_R_std_binorm_lcdf dist_fun;
  log_binorm tru_fun;
  Matrix<double, Dynamic, 1> inp_vec(3);
  inp_vec << 0.4, -2.7, 0.8;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun, inp_vec, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun, inp_vec, fx_fd, grad_fx_fd);

  EXPECT_FLOAT_EQ(fx_test, fx_tru);
  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i));
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_grad_test_vV_R) {
  int N_y = 3;
  vV_R_std_binorm_lcdf dist_fun2(N_y);
  vV_v_R_log_std_binorm_integral tru_fun(N_y);
  Matrix<double, Dynamic, 1> inp_vec2(N_y * 2 + 1);
  inp_vec2 << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7, 0.3;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun2, inp_vec2, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec2, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun2, inp_vec2, fx_fd, grad_fx_fd);

  EXPECT_FLOAT_EQ(fx_test, fx_tru);
  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i));
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_grad_test_vRV_R) {
  int N_y = 3;
  vRV_R_std_binorm_lcdf dist_fun2(N_y);
  vV_v_R_log_std_binorm_integral tru_fun(N_y);
  Matrix<double, Dynamic, 1> inp_vec2(N_y * 2 + 1);
  inp_vec2 << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7, 0.3;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun2, inp_vec2, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec2, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun2, inp_vec2, fx_fd, grad_fx_fd);

  EXPECT_FLOAT_EQ(fx_test, fx_tru);
  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i));
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_grad_test_vV_v) {
  int N_y = 3;
  vV_v_std_binorm_lcdf dist_fun(N_y);
  vV_v_R_log_std_binorm_integral tru_fun(N_y);
  Matrix<double, Dynamic, 1> inp_vec(N_y * 3);
  inp_vec << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7, 0.3, 0.4, 0.5;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun, inp_vec, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun, inp_vec, fx_fd, grad_fx_fd);

  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i));
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_grad_test_vRV_v) {
  int N_y = 3;
  vRV_v_std_binorm_lcdf dist_fun(N_y);
  vV_v_R_log_std_binorm_integral tru_fun(N_y);
  Matrix<double, Dynamic, 1> inp_vec(N_y * 3);
  inp_vec << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7, 0.3, 0.4, 0.5;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun, inp_vec, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun, inp_vec, fx_fd, grad_fx_fd);

  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i));
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_grad_test_vV_V) {
  int N_y = 3;
  vV_V_std_binorm_lcdf dist_fun(N_y);
  vV_v_R_log_std_binorm_integral tru_fun(N_y);
  Matrix<double, Dynamic, 1> inp_vec(N_y * 3);
  inp_vec << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7, 0.3, 0.4, 0.5;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun, inp_vec, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun, inp_vec, fx_fd, grad_fx_fd);

  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i));
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_grad_test_vRV_V) {
  int N_y = 3;
  vRV_V_std_binorm_lcdf dist_fun(N_y);
  vV_v_R_log_std_binorm_integral tru_fun(N_y);
  Matrix<double, Dynamic, 1> inp_vec(N_y * 3);
  inp_vec << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7, 0.3, 0.4, 0.5;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun, inp_vec, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun, inp_vec, fx_fd, grad_fx_fd);

  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i));
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_V_D_test) {
  V_D_std_binorm_lcdf dist_fun2(0.3);
  vector<double> rhos;
  rhos.push_back(0.3);
  vV_vD_D_log_std_binorm_integral tru_fun(1, rhos);
  Matrix<double, Dynamic, 1> inp_vec2(2);
  inp_vec2 << 0.4, -2.7;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun2, inp_vec2, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec2, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun2, inp_vec2, fx_fd, grad_fx_fd);

  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i)) << i;
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_RV_D_test) {
  RV_D_std_binorm_lcdf dist_fun2(0.3);
  vector<double> rhos;
  rhos.push_back(0.3);
  vV_vD_D_log_std_binorm_integral tru_fun(1, rhos);
  Matrix<double, Dynamic, 1> inp_vec2(2);
  inp_vec2 << 0.4, -2.7;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun2, inp_vec2, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec2, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun2, inp_vec2, fx_fd, grad_fx_fd);

  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i)) << i;
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_vV_D_test) {
  int N_y = 3;
  vV_D_std_binorm_lcdf dist_fun2(N_y, 0.3);
  vector<double> rhos;
  rhos.push_back(0.3);
  vV_vD_D_log_std_binorm_integral tru_fun(N_y, rhos);
  Matrix<double, Dynamic, 1> inp_vec2(N_y * 2);
  inp_vec2 << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun2, inp_vec2, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec2, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun2, inp_vec2, fx_fd, grad_fx_fd);

  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i)) << i;
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_vRV_D_test) {
  int N_y = 3;
  vRV_D_std_binorm_lcdf dist_fun2(N_y, 0.3);
  vector<double> rhos;
  rhos.push_back(0.3);
  vV_vD_D_log_std_binorm_integral tru_fun(N_y, rhos);
  Matrix<double, Dynamic, 1> inp_vec2(N_y * 2);
  inp_vec2 << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun2, inp_vec2, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec2, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun2, inp_vec2, fx_fd, grad_fx_fd);

  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i)) << i;
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_vV_vD_test) {
  int N_y = 3;
  vector<double> rhos;
  rhos.push_back(0.3);
  rhos.push_back(0.4);
  rhos.push_back(0.99);
  vV_vD_std_binorm_lcdf dist_fun2(N_y, rhos);
  vV_vD_D_log_std_binorm_integral tru_fun(N_y, rhos);
  Matrix<double, Dynamic, 1> inp_vec2(N_y * 2);
  inp_vec2 << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun2, inp_vec2, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec2, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun2, inp_vec2, fx_fd, grad_fx_fd);

  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i)) << i;
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_vRV_vD_test) {
  int N_y = 3;
  vector<double> rhos;
  rhos.push_back(0.3);
  rhos.push_back(0.4);
  rhos.push_back(0.99);
  vRV_vD_std_binorm_lcdf dist_fun2(N_y, rhos);
  vV_vD_D_log_std_binorm_integral tru_fun(N_y, rhos);
  Matrix<double, Dynamic, 1> inp_vec2(N_y * 2);
  inp_vec2 << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun2, inp_vec2, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec2, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun2, inp_vec2, fx_fd, grad_fx_fd);

  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i)) << i;
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_vV_VD_test) {
  int N_y = 3;
  Matrix<double, Dynamic, 1> rhos_eig(3);
  rhos_eig(0) = 0.3;
  rhos_eig(1) = 0.4;
  rhos_eig(2) = 0.99;
  vector<double> rhos;
  for (int i = 0; i < rhos_eig.size(); ++i)
    rhos.push_back(rhos_eig(i));
  vV_VD_std_binorm_lcdf dist_fun2(N_y, rhos_eig);
  vV_vD_D_log_std_binorm_integral tru_fun(N_y, rhos);
  Matrix<double, Dynamic, 1> inp_vec2(N_y * 2);
  inp_vec2 << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun2, inp_vec2, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec2, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun2, inp_vec2, fx_fd, grad_fx_fd);

  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i)) << i;
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_vRV_VD_test) {
  int N_y = 3;
  Matrix<double, Dynamic, 1> rhos_eig(3);
  rhos_eig(0) = 0.3;
  rhos_eig(1) = 0.4;
  rhos_eig(2) = 0.99;
  vector<double> rhos;
  for (int i = 0; i < rhos_eig.size(); ++i)
    rhos.push_back(rhos_eig(i));
  vRV_VD_std_binorm_lcdf dist_fun2(N_y, rhos_eig);
  vV_vD_D_log_std_binorm_integral tru_fun(N_y, rhos);
  Matrix<double, Dynamic, 1> inp_vec2(N_y * 2);
  inp_vec2 << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun2, inp_vec2, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec2, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun2, inp_vec2, fx_fd, grad_fx_fd);

  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i)) << i;
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_VD_R_test) {
  Matrix<double, Dynamic, 1> y_el(2);
  y_el << 2, 3;
  VD_R_std_binorm_lcdf dist_fun2(y_el);
  vector<Matrix<double, Dynamic, 1>> y;
  y.push_back(y_el);
  vVD_v_log_std_binorm_integral<Dynamic, 1> tru_fun(y);
  Matrix<double, Dynamic, 1> inp_vec2(1);
  inp_vec2(0) = 0.4;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun2, inp_vec2, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec2, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun2, inp_vec2, fx_fd, grad_fx_fd);

  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i)) << i;
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_RVD_R_test) {
  Matrix<double, 1, Dynamic> y_el(2);
  y_el << 2, 3;
  RVD_R_std_binorm_lcdf dist_fun2(y_el);
  vector<Matrix<double, 1, Dynamic>> y;
  y.push_back(y_el);
  vVD_v_log_std_binorm_integral<1, Dynamic> tru_fun(y);
  Matrix<double, Dynamic, 1> inp_vec2(1);
  inp_vec2(0) = 0.4;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun2, inp_vec2, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec2, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun2, inp_vec2, fx_fd, grad_fx_fd);

  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i)) << i;
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_vVD_R_test) {
  int N_y = 3;
  vector<Matrix<double, Dynamic, 1>> y;
  Matrix<double, Dynamic, 1> y_el(2);
  y_el << 2, 3;
  for (int i = 0; i < N_y; ++i)
    y.push_back(y_el);
  vVD_R_std_binorm_lcdf dist_fun2(y);
  vVD_v_log_std_binorm_integral<Dynamic, 1> tru_fun(y);
  Matrix<double, Dynamic, 1> inp_vec2(1);
  inp_vec2(0) = 0.4;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun2, inp_vec2, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec2, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun2, inp_vec2, fx_fd, grad_fx_fd);

  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i)) << i;
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_vRVD_R_test) {
  int N_y = 3;
  vector<Matrix<double, 1, Dynamic>> y;
  Matrix<double, 1, Dynamic> y_el(2);
  y_el << 2, 3;
  for (int i = 0; i < N_y; ++i)
    y.push_back(y_el);
  vRVD_R_std_binorm_lcdf dist_fun2(y);
  vVD_v_log_std_binorm_integral<1, Dynamic> tru_fun(y);
  Matrix<double, Dynamic, 1> inp_vec2(1);
  inp_vec2(0) = 0.4;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun2, inp_vec2, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec2, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun2, inp_vec2, fx_fd, grad_fx_fd);

  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i)) << i;
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_vVD_v_test) {
  int N_y = 3;
  vector<Matrix<double, Dynamic, 1>> y;
  Matrix<double, Dynamic, 1> y_el(2);
  y_el << 2, 3;
  for (int i = 0; i < N_y; ++i)
    y.push_back(y_el);
  vVD_v_std_binorm_lcdf dist_fun2(y);
  vVD_v_log_std_binorm_integral<Dynamic, 1> tru_fun(y);
  Matrix<double, Dynamic, 1> inp_vec2(3);
  inp_vec2 << 0.4, 0.3, -0.3;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun2, inp_vec2, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec2, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun2, inp_vec2, fx_fd, grad_fx_fd);

  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i)) << i;
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_vRVD_v_test) {
  int N_y = 3;
  vector<Matrix<double, 1, Dynamic>> y;
  Matrix<double, 1, Dynamic> y_el(2);
  y_el << 2, 3;
  for (int i = 0; i < N_y; ++i)
    y.push_back(y_el);
  vRVD_v_std_binorm_lcdf dist_fun2(y);
  vVD_v_log_std_binorm_integral<1, Dynamic> tru_fun(y);
  Matrix<double, Dynamic, 1> inp_vec2(3);
  inp_vec2 << 0.4, 0.3, -0.3;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun2, inp_vec2, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec2, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun2, inp_vec2, fx_fd, grad_fx_fd);

  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i)) << i;
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_vVD_V_test) {
  int N_y = 3;
  vector<Matrix<double, Dynamic, 1>> y;
  Matrix<double, Dynamic, 1> y_el(2);
  y_el << 2, 3;
  for (int i = 0; i < N_y; ++i)
    y.push_back(y_el);
  vVD_V_std_binorm_lcdf dist_fun2(y);
  vVD_v_log_std_binorm_integral<Dynamic, 1> tru_fun(y);
  Matrix<double, Dynamic, 1> inp_vec2(3);
  inp_vec2 << 0.4, 0.3, -0.3;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  stan::math::gradient(dist_fun2, inp_vec2, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec2, fx_tru, grad_fx_tru);
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::finite_diff_gradient(dist_fun2, inp_vec2, fx_fd, grad_fx_fd);

  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i)) << i;
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
TEST(MathFunctions, vec_binormal_integral_vRVD_V_test) {
  int N_y = 3;
  vector<Matrix<double, 1, Dynamic>> y;
  Matrix<double, 1, Dynamic> y_el(2);
  y_el << 2, 3;
  for (int i = 0; i < N_y; ++i)
    y.push_back(y_el);
  vRVD_V_std_binorm_lcdf dist_fun2(y);
  vVD_v_log_std_binorm_integral<1, Dynamic> tru_fun(y);
  Matrix<double, Dynamic, 1> inp_vec2(3);
  inp_vec2 << 0.4, 0.3, -0.3;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  double fx_fd;
  Matrix<double, Dynamic, 1> grad_fx_fd;
  stan::math::gradient(dist_fun2, inp_vec2, fx_test, grad_fx_test);
  stan::math::gradient(tru_fun, inp_vec2, fx_tru, grad_fx_tru);
  stan::math::finite_diff_gradient(dist_fun2, inp_vec2, fx_fd, grad_fx_fd);

  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i));
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_fd(i));
  }
}
