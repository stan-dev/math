#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <limits>
#include <vector>
#include <tuple>

using Eigen::Dynamic;
using Eigen::Matrix;
using stan::math::var;
using std::vector;

// Vector/Row Vector, real
template <int R, int C>
struct V_R_binorm_copula_cdf {
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    Matrix<T, R, C> y(2);
    y(0) = inp_vec(0);
    y(1) = inp_vec(1);
    return stan::math::binormal_copula_cdf(y, inp_vec(2));
  }
};

// Vector/Row Vector, double
template <int R, int C>
struct V_D_binorm_copula_cdf {
  double rho_;
  explicit V_D_binorm_copula_cdf(double rho) : rho_(rho) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    Matrix<T, R, C> y(2);
    y(0) = inp_vec(0);
    y(1) = inp_vec(1);
    return stan::math::binormal_copula_cdf(y, rho_);
  }
};

// double Vector/Row Vector, real
template <int R, int C>
struct VD_R_binorm_copula_cdf {
  Matrix<double, R, C> y_;
  explicit VD_R_binorm_copula_cdf(const Matrix<double, R, C>& y) : y_(y) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    return stan::math::binormal_copula_cdf(y_, inp_vec(0));
  }
};


struct V_R_binorm_copula_manual {
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    using std::log;
    using stan::math::inv_Phi;
    return
        stan::math::std_binormal_integral(inv_Phi(inp_vec(0)),
                                          inv_Phi(inp_vec(1)), inp_vec(2));
  }
};

struct V_D_binorm_copula_manual {
  double rho_;
  explicit V_D_binorm_copula_manual(const double rho) : rho_(rho) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    using std::log;
    using stan::math::inv_Phi;
    return
        stan::math::std_binormal_integral(inv_Phi(inp_vec(0)),
                                          inv_Phi(inp_vec(1)), rho_);
  }
};

template <int R, int C>
struct VD_R_binorm_copula_manual {
  Matrix<double, R, C> y_;
  explicit VD_R_binorm_copula_manual(const Matrix<double, R, C>& y) : y_(y) {}
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    using std::log;
    using stan::math::inv_Phi;
    return
        stan::math::std_binormal_integral(inv_Phi(y_(0)),
                                          inv_Phi(y_(1)), inp_vec(0));
  }
};

template<typename T1, typename T2, int R, int C>
typename stan::return_type<T1, T2>::type
call_args(const std::tuple<Matrix<T1, R, C>, T2>& inps) {
  auto arg1 = std::get<0>(inps);
  auto arg2 = std::get<1>(inps);
  return stan::math::binormal_copula_cdf(arg1, arg2);
}

template<typename T1, typename T2, int R, int C>
std::tuple<Matrix<T1, R, C>, T2>
make_args(const Matrix<double, Dynamic, 1>& flat_args) {
  int dim_y = flat_args.size() - 1;
  Eigen::Matrix<T1, R, C> y(dim_y);
  for (int i = 0; i < dim_y; ++i) {
    y(i) = T1(flat_args(i));
  }

  T2 rho = T2(flat_args(dim_y));
  return std::make_tuple(y, rho);
}
