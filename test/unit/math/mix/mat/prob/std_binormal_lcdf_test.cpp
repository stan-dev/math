#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <limits>
#include <vector>

using Eigen::Dynamic;
using Eigen::Matrix;
using stan::math::fvar;
using stan::math::var;
using std::vector;

// Vector, real
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

TEST(MathFunctions, binormal_integral_throw_RV_1_nan_fv) {
  var nan = std::numeric_limits<double>::quiet_NaN();
  Matrix<fvar<var>, Dynamic, 1> y(2);
  y(0) = fvar<var>(var(nan), var(0));
  y(1) = fvar<var>(var(-2.0), var(0));
  fvar<var> rho(0.3, 0);
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_RV_2_nan_fv) {
  var nan = std::numeric_limits<double>::quiet_NaN();
  Matrix<fvar<var>, Dynamic, 1> y(2);
  y(0) = fvar<var>(var(-2.0), var(0));
  y(1) = fvar<var>(var(nan), var(0));
  fvar<var> rho(0.3, 0);
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_nan_fv) {
  var nan = std::numeric_limits<double>::quiet_NaN();
  Matrix<fvar<var>, Dynamic, 1> y(2);
  y(0) = fvar<var>(var(-2.0), var(0));
  y(1) = fvar<var>(var(1), var(0));
  fvar<var> rho(nan, 0);
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_LT_neg_one_fv) {
  Matrix<fvar<var>, Dynamic, 1> y(2);
  y(0) = fvar<var>(var(-2.0), var(0));
  y(1) = fvar<var>(var(1), var(0));
  fvar<var> rho(-1.3, 0);
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_GT_one_fv) {
  Matrix<fvar<var>, Dynamic, 1> y(2);
  y(0) = fvar<var>(var(-2.0), var(0));
  y(1) = fvar<var>(var(1), var(0));
  fvar<var> rho(1.3, 0);
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_RV_1_nan_ffv) {
  fvar<var> nan(std::numeric_limits<double>::quiet_NaN(), 0);
  Matrix<fvar<fvar<var>>, Dynamic, 1> y(2);
  y(0) = fvar<fvar<var>>(nan, fvar<var>(3, 0));
  y(1) = fvar<fvar<var>>(fvar<var>(2, 0), fvar<var>(3, 0));
  fvar<fvar<var>> rho(fvar<var>(0.3, 0), fvar<var>(0, 0));
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_RV_2_nan_ffv) {
  fvar<var> nan(std::numeric_limits<double>::quiet_NaN(), 0);
  Matrix<fvar<fvar<var>>, Dynamic, 1> y(2);
  y(0) = fvar<fvar<var>>(fvar<var>(2, 0), fvar<var>(3, 0));
  y(1) = fvar<fvar<var>>(nan, fvar<var>(3, 0));
  fvar<fvar<var>> rho(fvar<var>(0.3, 0), fvar<var>(0, 0));
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_nan_ffv) {
  fvar<var> nan(std::numeric_limits<double>::quiet_NaN(), 0);
  Matrix<fvar<fvar<var>>, Dynamic, 1> y(2);
  y(0) = fvar<fvar<var>>(fvar<var>(2, 0), fvar<var>(3, 0));
  y(1) = fvar<fvar<var>>(fvar<var>(3, 0), fvar<var>(3, 0));
  fvar<fvar<var>> rho(nan, fvar<var>(0, 0));
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_LT_neg_one_ffv) {
  Matrix<fvar<fvar<var>>, Dynamic, 1> y(2);
  y(0) = fvar<fvar<var>>(fvar<var>(2, 0), fvar<var>(3, 0));
  y(1) = fvar<fvar<var>>(fvar<var>(3, 0), fvar<var>(3, 0));
  fvar<fvar<var>> rho(fvar<var>(-1.3, 0), fvar<var>(0, 0));
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_GT_one_ffv) {
  Matrix<fvar<fvar<var>>, Dynamic, 1> y(2);
  y(0) = fvar<fvar<var>>(fvar<var>(2, 0), fvar<var>(3, 0));
  y(1) = fvar<fvar<var>>(fvar<var>(3, 0), fvar<var>(3, 0));
  fvar<fvar<var>> rho(fvar<var>(1.3, 0), fvar<var>(0, 0));
  EXPECT_THROW(stan::math::std_binormal_lcdf(y, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_val_test_fv) {
  vector<Matrix<fvar<var>, Dynamic, 1>> vals(14,
                                             Matrix<fvar<var>, Dynamic, 1>(3));
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
    EXPECT_FLOAT_EQ(dist_fun(vals[i]).val_.val(), tru_fun(vals[i]).val_.val());
}
TEST(MathFunctions, binormal_integral_val_test_ffv) {
  vector<Matrix<fvar<fvar<var>>, Dynamic, 1>> vals(
      14, Matrix<fvar<fvar<var>>, Dynamic, 1>(3));
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
    EXPECT_FLOAT_EQ(dist_fun(vals[i]).val_.val_.val(),
                    tru_fun(vals[i]).val_.val_.val());
}
TEST(MathFunctions, vec_binormal_integral_val_test_vV_v_fv) {
  int N_y = 3;
  vV_v_std_binorm_lcdf dist_fun(N_y);
  vV_R_std_binorm_lcdf dist_fun2(N_y);
  vV_v_R_log_std_binorm_integral tru_fun(N_y);
  Matrix<fvar<var>, Dynamic, 1> inp_vec(N_y * 3);
  inp_vec << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7, 0.3, 0.4, 0.5;

  EXPECT_FLOAT_EQ(dist_fun(inp_vec).val_.val(), tru_fun(inp_vec).val_.val());

  Matrix<fvar<var>, Dynamic, 1> inp_vec2(N_y * 2 + 1);
  inp_vec2 << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7, 0.3;
  EXPECT_FLOAT_EQ(dist_fun2(inp_vec2).val_.val(), tru_fun(inp_vec2).val_.val());
}
TEST(MathFunctions, vec_binormal_integral_val_test_vV_v_ffv) {
  int N_y = 3;
  vV_v_std_binorm_lcdf dist_fun(N_y);
  vV_R_std_binorm_lcdf dist_fun2(N_y);
  vV_v_R_log_std_binorm_integral tru_fun(N_y);
  Matrix<fvar<fvar<var>>, Dynamic, 1> inp_vec(N_y * 3);
  inp_vec << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7, 0.3, 0.4, 0.5;

  EXPECT_FLOAT_EQ(dist_fun(inp_vec).val_.val_.val(),
                  tru_fun(inp_vec).val_.val_.val());

  Matrix<fvar<fvar<var>>, Dynamic, 1> inp_vec2(N_y * 2 + 1);
  inp_vec2 << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7, 0.3;
  EXPECT_FLOAT_EQ(dist_fun2(inp_vec2).val_.val_.val(),
                  tru_fun(inp_vec2).val_.val_.val());
}
TEST(MathFunctions, vec_binormal_integral_hess_test_V_R) {
  V_R_std_binorm_lcdf dist_fun;
  log_binorm tru_fun;
  Matrix<double, Dynamic, 1> inp_vec(3);
  inp_vec << 0.4, -2.7, 0.9;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  Matrix<double, Dynamic, Dynamic> H_test;

  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  Matrix<double, Dynamic, Dynamic> H_tru;
  stan::math::hessian(dist_fun, inp_vec, fx_test, grad_fx_test, H_test);
  stan::math::hessian(tru_fun, inp_vec, fx_tru, grad_fx_tru, H_tru);

  EXPECT_FLOAT_EQ(fx_test, fx_tru);
  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i));
  }
  for (int i = 0; i < grad_fx_test.size(); ++i)
    for (int j = 0; j < grad_fx_test.size(); ++j)
      EXPECT_FLOAT_EQ(H_test(i, j), H_tru(i, j));
}
TEST(MathFunctions, vec_binormal_integral_hess_test_vV_R) {
  int N_y = 3;
  vV_R_std_binorm_lcdf dist_fun2(N_y);
  vV_v_R_log_std_binorm_integral tru_fun(N_y);
  Matrix<double, Dynamic, 1> inp_vec2(N_y * 2 + 1);
  inp_vec2 << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7, 0.3;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  Matrix<double, Dynamic, Dynamic> H_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  Matrix<double, Dynamic, Dynamic> H_tru;
  stan::math::hessian(dist_fun2, inp_vec2, fx_test, grad_fx_test, H_test);
  stan::math::hessian(tru_fun, inp_vec2, fx_tru, grad_fx_tru, H_tru);

  EXPECT_FLOAT_EQ(fx_test, fx_tru);
  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i));
  }
  for (int i = 0; i < H_tru.rows(); ++i)
    for (int j = 0; j < H_tru.cols(); ++j)
      EXPECT_FLOAT_EQ(H_test(i, j), H_tru(i, j));
}
TEST(MathFunctions, vec_binormal_integral_hess_test_vV_v) {
  int N_y = 3;
  vV_v_std_binorm_lcdf dist_fun2(N_y);
  vV_v_R_log_std_binorm_integral tru_fun(N_y);
  Matrix<double, Dynamic, 1> inp_vec(N_y * 3);
  inp_vec << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7, 0.3, 0.4, 0.5;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  Matrix<double, Dynamic, Dynamic> H_test;
  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  Matrix<double, Dynamic, Dynamic> H_tru;
  stan::math::hessian(dist_fun2, inp_vec, fx_test, grad_fx_test, H_test);
  stan::math::hessian(tru_fun, inp_vec, fx_tru, grad_fx_tru, H_tru);

  EXPECT_FLOAT_EQ(fx_test, fx_tru);
  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i));
  }
  for (int i = 0; i < H_tru.rows(); ++i)
    for (int j = 0; j < H_tru.cols(); ++j)
      EXPECT_FLOAT_EQ(H_test(i, j), H_tru(i, j));
}
TEST(MathFunctions, vec_binormal_integral_hess_test_V_D) {
  V_D_std_binorm_lcdf dist_fun(0.3);
  vector<double> rho;
  rho.push_back(0.3);
  vV_vD_D_log_std_binorm_integral tru_fun(1, rho);
  Matrix<double, Dynamic, 1> inp_vec(2);
  inp_vec << 0.4, -2.7;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  Matrix<double, Dynamic, Dynamic> H_test;

  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  Matrix<double, Dynamic, Dynamic> H_tru;
  stan::math::hessian(dist_fun, inp_vec, fx_test, grad_fx_test, H_test);
  stan::math::hessian(tru_fun, inp_vec, fx_tru, grad_fx_tru, H_tru);

  EXPECT_FLOAT_EQ(fx_test, fx_tru);
  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i));
  }
  for (int i = 0; i < grad_fx_test.size(); ++i)
    for (int j = 0; j < grad_fx_test.size(); ++j)
      EXPECT_FLOAT_EQ(H_test(i, j), H_tru(i, j));
}
TEST(MathFunctions, vec_binormal_integral_hess_test_vV_D) {
  vV_D_std_binorm_lcdf dist_fun(3, 0.3);
  vector<double> rho;
  rho.push_back(0.3);
  vV_vD_D_log_std_binorm_integral tru_fun(3, rho);
  Matrix<double, Dynamic, 1> inp_vec(6);
  inp_vec << 0.4, -2.7, 0.4, 2, 2, -1.5;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  Matrix<double, Dynamic, Dynamic> H_test;

  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  Matrix<double, Dynamic, Dynamic> H_tru;
  stan::math::hessian(dist_fun, inp_vec, fx_test, grad_fx_test, H_test);
  stan::math::hessian(tru_fun, inp_vec, fx_tru, grad_fx_tru, H_tru);

  EXPECT_FLOAT_EQ(fx_test, fx_tru);
  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i));
  }
  for (int i = 0; i < grad_fx_test.size(); ++i)
    for (int j = 0; j < grad_fx_test.size(); ++j)
      EXPECT_FLOAT_EQ(H_test(i, j), H_tru(i, j));
}
TEST(MathFunctions, vec_binormal_integral_hess_test_vV_vD) {
  vector<double> rho;
  rho.push_back(0.3);
  rho.push_back(-0.3);
  rho.push_back(0.9);
  vV_vD_std_binorm_lcdf dist_fun(3, rho);
  vV_vD_D_log_std_binorm_integral tru_fun(3, rho);
  Matrix<double, Dynamic, 1> inp_vec(6);
  inp_vec << 0.4, -2.7, 0.4, 2, 2, -1.5;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  Matrix<double, Dynamic, Dynamic> H_test;

  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  Matrix<double, Dynamic, Dynamic> H_tru;
  stan::math::hessian(dist_fun, inp_vec, fx_test, grad_fx_test, H_test);
  stan::math::hessian(tru_fun, inp_vec, fx_tru, grad_fx_tru, H_tru);

  EXPECT_FLOAT_EQ(fx_test, fx_tru);
  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i));
  }
  for (int i = 0; i < grad_fx_test.size(); ++i)
    for (int j = 0; j < grad_fx_test.size(); ++j)
      EXPECT_FLOAT_EQ(H_test(i, j), H_tru(i, j));
}
TEST(MathFunctions, vec_binormal_integral_hess_test_VD_R) {
  Matrix<double, Dynamic, 1> y_el(2);
  y_el << 2, 3;
  vector<Matrix<double, Dynamic, 1>> y;
  y.push_back(y_el);
  VD_R_std_binorm_lcdf dist_fun(y_el);
  vVD_v_log_std_binorm_integral<Dynamic, 1> tru_fun(y);
  Matrix<double, Dynamic, 1> inp_vec(1);
  inp_vec << 0.4;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  Matrix<double, Dynamic, Dynamic> H_test;

  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  Matrix<double, Dynamic, Dynamic> H_tru;
  stan::math::hessian(dist_fun, inp_vec, fx_test, grad_fx_test, H_test);
  stan::math::hessian(tru_fun, inp_vec, fx_tru, grad_fx_tru, H_tru);

  EXPECT_FLOAT_EQ(fx_test, fx_tru);
  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i));
  }
  for (int i = 0; i < grad_fx_test.size(); ++i)
    for (int j = 0; j < grad_fx_test.size(); ++j)
      EXPECT_FLOAT_EQ(H_test(i, j), H_tru(i, j));
}
TEST(MathFunctions, vec_binormal_integral_hess_test_vVD_R) {
  Matrix<double, Dynamic, 1> y_el(2);
  y_el << 2, 3;
  vector<Matrix<double, Dynamic, 1>> y;
  y.push_back(y_el);
  y.push_back(y_el);
  y.push_back(y_el);
  vVD_R_std_binorm_lcdf dist_fun(y);
  vVD_v_log_std_binorm_integral<Dynamic, 1> tru_fun(y);
  Matrix<double, Dynamic, 1> inp_vec(1);
  inp_vec << 0.4;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  Matrix<double, Dynamic, Dynamic> H_test;

  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  Matrix<double, Dynamic, Dynamic> H_tru;
  stan::math::hessian(dist_fun, inp_vec, fx_test, grad_fx_test, H_test);
  stan::math::hessian(tru_fun, inp_vec, fx_tru, grad_fx_tru, H_tru);

  EXPECT_FLOAT_EQ(fx_test, fx_tru);
  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i));
  }
  for (int i = 0; i < grad_fx_test.size(); ++i)
    for (int j = 0; j < grad_fx_test.size(); ++j)
      EXPECT_FLOAT_EQ(H_test(i, j), H_tru(i, j));
}
TEST(MathFunctions, vec_binormal_integral_hess_test_vVD_v) {
  Matrix<double, Dynamic, 1> y_el(2);
  y_el << 2, 3;
  vector<Matrix<double, Dynamic, 1>> y;
  y.push_back(y_el);
  y.push_back(y_el);
  y.push_back(y_el);
  vVD_v_std_binorm_lcdf dist_fun(y);
  vVD_v_log_std_binorm_integral<Dynamic, 1> tru_fun(y);
  Matrix<double, Dynamic, 1> inp_vec(3);
  inp_vec << 0.4, -0.3, 0.9;

  double fx_test;
  Matrix<double, Dynamic, 1> grad_fx_test;
  Matrix<double, Dynamic, Dynamic> H_test;

  double fx_tru;
  Matrix<double, Dynamic, 1> grad_fx_tru;
  Matrix<double, Dynamic, Dynamic> H_tru;
  stan::math::hessian(dist_fun, inp_vec, fx_test, grad_fx_test, H_test);
  stan::math::hessian(tru_fun, inp_vec, fx_tru, grad_fx_tru, H_tru);

  EXPECT_FLOAT_EQ(fx_test, fx_tru);
  for (int i = 0; i < grad_fx_test.size(); ++i) {
    EXPECT_FLOAT_EQ(grad_fx_test(i), grad_fx_tru(i));
  }
  for (int i = 0; i < grad_fx_test.size(); ++i)
    for (int j = 0; j < grad_fx_test.size(); ++j)
      EXPECT_FLOAT_EQ(H_test(i, j), H_tru(i, j));
}
TEST(MathFunctions, vec_binormal_integral_grad_hess_test_V_R) {
  V_R_std_binorm_lcdf dist_fun;
  log_binorm tru_fun;
  Matrix<double, Dynamic, 1> inp_vec(3);
  inp_vec << 0.4, -2.7, 0.9;

  double fx_test;
  Matrix<double, Dynamic, Dynamic> H_test;

  double fx_tru;
  Matrix<double, Dynamic, Dynamic> H_tru;
  vector<Matrix<double, Dynamic, Dynamic>> grad_H_test;
  vector<Matrix<double, Dynamic, Dynamic>> grad_H_tru;
  stan::math::grad_hessian(dist_fun, inp_vec, fx_test, H_test, grad_H_test);
  stan::math::grad_hessian(tru_fun, inp_vec, fx_tru, H_tru, grad_H_tru);

  EXPECT_FLOAT_EQ(fx_test, fx_tru);
  for (size_t i = 0; i < grad_H_test.size(); ++i)
    for (int j = 0; j < grad_H_test[i].rows(); ++j)
      for (int k = 0; k < grad_H_test[i].cols(); ++k)
        EXPECT_FLOAT_EQ(grad_H_test[i](j, k), grad_H_tru[i](j, k));
}
TEST(MathFunctions, vec_binormal_integral_grad_hess_test_vV_R) {
  int N_y = 3;
  vV_R_std_binorm_lcdf dist_fun2(N_y);
  vV_v_R_log_std_binorm_integral tru_fun(N_y);
  Matrix<double, Dynamic, 1> inp_vec2(N_y * 2 + 1);
  inp_vec2 << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7, 0.3;

  double fx_test;
  Matrix<double, Dynamic, Dynamic> H_test;
  double fx_tru;
  Matrix<double, Dynamic, Dynamic> H_tru;
  vector<Matrix<double, Dynamic, Dynamic>> grad_H_test;
  vector<Matrix<double, Dynamic, Dynamic>> grad_H_tru;
  stan::math::grad_hessian(dist_fun2, inp_vec2, fx_test, H_test, grad_H_test);
  stan::math::grad_hessian(tru_fun, inp_vec2, fx_tru, H_tru, grad_H_tru);

  EXPECT_FLOAT_EQ(fx_test, fx_tru);
  for (size_t i = 0; i < grad_H_test.size(); ++i)
    for (int j = 0; j < grad_H_test[i].rows(); ++j)
      for (int k = 0; k < grad_H_test[i].cols(); ++k)
        EXPECT_FLOAT_EQ(grad_H_test[i](j, k), grad_H_tru[i](j, k));
}
TEST(MathFunctions, vec_binormal_integral_grad_hess_test_vV_v) {
  int N_y = 3;
  vV_v_std_binorm_lcdf dist_fun2(N_y);
  vV_v_R_log_std_binorm_integral tru_fun(N_y);
  Matrix<double, Dynamic, 1> inp_vec(N_y * 3);
  inp_vec << 0.4, -2.7, 0.4, -2.7, 0.4, -2.7, 0.3, 0.4, 0.5;

  double fx_test;
  Matrix<double, Dynamic, Dynamic> H_test;
  double fx_tru;
  Matrix<double, Dynamic, Dynamic> H_tru;
  vector<Matrix<double, Dynamic, Dynamic>> grad_H_test;
  vector<Matrix<double, Dynamic, Dynamic>> grad_H_tru;
  stan::math::grad_hessian(dist_fun2, inp_vec, fx_test, H_test, grad_H_test);
  stan::math::grad_hessian(tru_fun, inp_vec, fx_tru, H_tru, grad_H_tru);

  EXPECT_FLOAT_EQ(fx_test, fx_tru);
  for (size_t i = 0; i < grad_H_test.size(); ++i)
    for (int j = 0; j < grad_H_test[i].rows(); ++j)
      for (int k = 0; k < grad_H_test[i].cols(); ++k)
        EXPECT_FLOAT_EQ(grad_H_test[i](j, k), grad_H_tru[i](j, k));
}
TEST(MathFunctions, vec_binormal_integral_grad_hess_test_V_D) {
  V_D_std_binorm_lcdf dist_fun(0.3);
  vector<double> rho;
  rho.push_back(0.3);
  vV_vD_D_log_std_binorm_integral tru_fun(1, rho);
  Matrix<double, Dynamic, 1> inp_vec(2);
  inp_vec << 0.4, -2.7;

  double fx_test;
  Matrix<double, Dynamic, Dynamic> H_test;

  double fx_tru;
  Matrix<double, Dynamic, Dynamic> H_tru;
  vector<Matrix<double, Dynamic, Dynamic>> grad_H_test;
  vector<Matrix<double, Dynamic, Dynamic>> grad_H_tru;
  stan::math::grad_hessian(dist_fun, inp_vec, fx_test, H_test, grad_H_test);
  stan::math::grad_hessian(tru_fun, inp_vec, fx_tru, H_tru, grad_H_tru);

  EXPECT_FLOAT_EQ(fx_test, fx_tru);
  for (size_t i = 0; i < grad_H_test.size(); ++i)
    for (int j = 0; j < grad_H_test[i].rows(); ++j)
      for (int k = 0; k < grad_H_test[i].cols(); ++k)
        EXPECT_FLOAT_EQ(grad_H_test[i](j, k), grad_H_tru[i](j, k));
}
TEST(MathFunctions, vec_binormal_integral_grad_hess_test_vV_D) {
  vV_D_std_binorm_lcdf dist_fun(3, 0.3);
  vector<double> rho;
  rho.push_back(0.3);
  vV_vD_D_log_std_binorm_integral tru_fun(3, rho);
  Matrix<double, Dynamic, 1> inp_vec(6);
  inp_vec << 0.4, -2.7, 0.4, 2, 2, -1.5;

  double fx_test;
  Matrix<double, Dynamic, Dynamic> H_test;

  double fx_tru;
  Matrix<double, Dynamic, Dynamic> H_tru;
  vector<Matrix<double, Dynamic, Dynamic>> grad_H_test;
  vector<Matrix<double, Dynamic, Dynamic>> grad_H_tru;
  stan::math::grad_hessian(dist_fun, inp_vec, fx_test, H_test, grad_H_test);
  stan::math::grad_hessian(tru_fun, inp_vec, fx_tru, H_tru, grad_H_tru);

  EXPECT_FLOAT_EQ(fx_test, fx_tru);
  for (size_t i = 0; i < grad_H_test.size(); ++i)
    for (int j = 0; j < grad_H_test[i].rows(); ++j)
      for (int k = 0; k < grad_H_test[i].cols(); ++k)
        EXPECT_FLOAT_EQ(grad_H_test[i](j, k), grad_H_tru[i](j, k));
}
TEST(MathFunctions, vec_binormal_integral_grad_hess_test_vV_vD) {
  vector<double> rho;
  rho.push_back(0.3);
  rho.push_back(-0.3);
  rho.push_back(0.9);
  vV_vD_std_binorm_lcdf dist_fun(3, rho);
  vV_vD_D_log_std_binorm_integral tru_fun(3, rho);
  Matrix<double, Dynamic, 1> inp_vec(6);
  inp_vec << 0.4, -2.7, 0.4, 2, 2, -1.5;

  double fx_test;
  Matrix<double, Dynamic, Dynamic> H_test;

  double fx_tru;
  Matrix<double, Dynamic, Dynamic> H_tru;
  vector<Matrix<double, Dynamic, Dynamic>> grad_H_test;
  vector<Matrix<double, Dynamic, Dynamic>> grad_H_tru;
  stan::math::grad_hessian(dist_fun, inp_vec, fx_test, H_test, grad_H_test);
  stan::math::grad_hessian(tru_fun, inp_vec, fx_tru, H_tru, grad_H_tru);

  EXPECT_FLOAT_EQ(fx_test, fx_tru);
  for (size_t i = 0; i < grad_H_test.size(); ++i)
    for (int j = 0; j < grad_H_test[i].rows(); ++j)
      for (int k = 0; k < grad_H_test[i].cols(); ++k)
        EXPECT_FLOAT_EQ(grad_H_test[i](j, k), grad_H_tru[i](j, k));
}
TEST(MathFunctions, vec_binormal_integral_grad_hess_test_VD_R) {
  Matrix<double, Dynamic, 1> y_el(2);
  y_el << 2, 3;
  vector<Matrix<double, Dynamic, 1>> y;
  y.push_back(y_el);
  VD_R_std_binorm_lcdf dist_fun(y_el);
  vVD_v_log_std_binorm_integral<Dynamic, 1> tru_fun(y);
  Matrix<double, Dynamic, 1> inp_vec(1);
  inp_vec << 0.4;

  double fx_test;
  Matrix<double, Dynamic, Dynamic> H_test;

  double fx_tru;
  Matrix<double, Dynamic, Dynamic> H_tru;
  vector<Matrix<double, Dynamic, Dynamic>> grad_H_test;
  vector<Matrix<double, Dynamic, Dynamic>> grad_H_tru;
  stan::math::grad_hessian(dist_fun, inp_vec, fx_test, H_test, grad_H_test);
  stan::math::grad_hessian(tru_fun, inp_vec, fx_tru, H_tru, grad_H_tru);

  EXPECT_FLOAT_EQ(fx_test, fx_tru);
  for (size_t i = 0; i < grad_H_test.size(); ++i)
    for (int j = 0; j < grad_H_test[i].rows(); ++j)
      for (int k = 0; k < grad_H_test[i].cols(); ++k)
        EXPECT_FLOAT_EQ(grad_H_test[i](j, k), grad_H_tru[i](j, k));
}
TEST(MathFunctions, vec_binormal_integral_grad_hess_test_vVD_R) {
  Matrix<double, Dynamic, 1> y_el(2);
  y_el << 2, 3;
  vector<Matrix<double, Dynamic, 1>> y;
  y.push_back(y_el);
  y.push_back(y_el);
  y.push_back(y_el);
  vVD_R_std_binorm_lcdf dist_fun(y);
  vVD_v_log_std_binorm_integral<Dynamic, 1> tru_fun(y);
  Matrix<double, Dynamic, 1> inp_vec(1);
  inp_vec << 0.4;

  double fx_test;
  Matrix<double, Dynamic, Dynamic> H_test;

  double fx_tru;
  Matrix<double, Dynamic, Dynamic> H_tru;
  vector<Matrix<double, Dynamic, Dynamic>> grad_H_test;
  vector<Matrix<double, Dynamic, Dynamic>> grad_H_tru;
  stan::math::grad_hessian(dist_fun, inp_vec, fx_test, H_test, grad_H_test);
  stan::math::grad_hessian(tru_fun, inp_vec, fx_tru, H_tru, grad_H_tru);

  EXPECT_FLOAT_EQ(fx_test, fx_tru);
  for (size_t i = 0; i < grad_H_test.size(); ++i)
    for (int j = 0; j < grad_H_test[i].rows(); ++j)
      for (int k = 0; k < grad_H_test[i].cols(); ++k)
        EXPECT_FLOAT_EQ(grad_H_test[i](j, k), grad_H_tru[i](j, k));
}
TEST(MathFunctions, vec_binormal_integral_grad_hess_test_vVD_v) {
  Matrix<double, Dynamic, 1> y_el(2);
  y_el << 2, 3;
  vector<Matrix<double, Dynamic, 1>> y;
  y.push_back(y_el);
  y.push_back(y_el);
  y.push_back(y_el);
  vVD_v_std_binorm_lcdf dist_fun(y);
  vVD_v_log_std_binorm_integral<Dynamic, 1> tru_fun(y);
  Matrix<double, Dynamic, 1> inp_vec(3);
  inp_vec << 0.4, -0.3, 0.9;

  double fx_test;
  Matrix<double, Dynamic, Dynamic> H_test;

  double fx_tru;
  Matrix<double, Dynamic, Dynamic> H_tru;
  vector<Matrix<double, Dynamic, Dynamic>> grad_H_test;
  vector<Matrix<double, Dynamic, Dynamic>> grad_H_tru;
  stan::math::grad_hessian(dist_fun, inp_vec, fx_test, H_test, grad_H_test);
  stan::math::grad_hessian(tru_fun, inp_vec, fx_tru, H_tru, grad_H_tru);

  EXPECT_FLOAT_EQ(fx_test, fx_tru);
  for (size_t i = 0; i < grad_H_test.size(); ++i)
    for (int j = 0; j < grad_H_test[i].rows(); ++j)
      for (int k = 0; k < grad_H_test[i].cols(); ++k)
        EXPECT_FLOAT_EQ(grad_H_test[i](j, k), grad_H_tru[i](j, k));
}
