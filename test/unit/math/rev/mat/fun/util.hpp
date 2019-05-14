#ifndef TEST_UNIT_MATH_REV_MAT_FUN_UTIL_HPP
#define TEST_UNIT_MATH_REV_MAT_FUN_UTIL_HPP

#include <test/unit/math/rev/arr/fun/util.hpp>
#include <stan/math/rev/mat.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <vector>

typedef stan::math::index_type<Eigen::Matrix<double, -1, -1> >::type size_type;

// Returns a matrix with the contents of a
// vector; Fills the matrix column-wise

template <typename T, int R, int C>
void fill(const std::vector<double>& contents, Eigen::Matrix<T, R, C>& M) {
  size_t ij = 0;
  for (int j = 0; j < C; ++j)
    for (int i = 0; i < R; ++i)
      M(i, j) = T(contents[ij++]);
}

template <typename T>
void create_vec(const std::vector<double>& vals, std::vector<T>& created_vars) {
  for (size_t i = 0; i < vals.size(); ++i)
    created_vars.push_back(T(vals[i]));
}

// fun3: R^3 --> R | (x, y, z) -- > x^3 * y^2 + x * y^3 + z^3 * x * y
struct third_order_mixed {
  template <typename T>
  inline T operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const {
    return x(0) * x(0) * x(0) * x(1) * x(1) + x(1) * x(1) * x(1) * x(0)
           + x(2) * x(2) * x(2) * x(1) * x(0);
  }
};

Eigen::Matrix<double, 3, 3> third_order_mixed_hess(
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& inp_vec) {
  Eigen::Matrix<double, 3, 3> hess;

  double x = inp_vec(0);
  double y = inp_vec(1);
  double z = inp_vec(2);

  double z_sq = z * z;
  double y_sq = y * y;
  double x_sq = x * x;
  double z_cub = z_sq * z;

  double f_xy = 6 * x_sq * y + 3 * y_sq + z_cub;

  hess << 6 * x * y_sq, f_xy, 3 * z_sq * y, f_xy, 2 * x_sq * x + 6 * x * y,
      3 * z_sq * x, 3 * z_sq * y, 3 * z_sq * x, 6 * x * y * z;
  return hess;
}

Eigen::Matrix<double, 3, 3> norm_hess(
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& inp_vec) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  Matrix<double, 3, 3> hess;
  double inv_sigma_sq = 1 / (inp_vec(2) * inp_vec(2));
  double y_m_mu = inp_vec(0) - inp_vec(1);
  double part_1_3 = 2 * y_m_mu * inv_sigma_sq / inp_vec(2);
  double part_3_3
      = inv_sigma_sq - 3 * inv_sigma_sq * inv_sigma_sq * y_m_mu * y_m_mu;
  hess << -inv_sigma_sq, inv_sigma_sq, part_1_3, inv_sigma_sq, -inv_sigma_sq,
      -part_1_3, part_1_3, -part_1_3, part_3_3;
  return hess;
}

std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >
third_order_mixed_grad_hess(
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& inp_vec) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  std::vector<Matrix<double, Dynamic, Dynamic> > grad_hess_ret;
  for (int i = 0; i < inp_vec.size(); ++i)
    grad_hess_ret.push_back(Matrix<double, Dynamic, Dynamic>(3, 3));

  double x = inp_vec(0);
  double y = inp_vec(1);
  double z = inp_vec(2);
  double x_sq = x * x;
  double y_sq = y * y;
  double z_sq = z * z;
  double zy = z * y;
  double zx = z * x;
  double yx = x * y;
  double xy = yx;

  grad_hess_ret[0] << 6 * y_sq, 12 * xy, 0, 12 * xy, 6 * x_sq + 6 * y, 3 * z_sq,
      0, 3 * z_sq, 6 * zy;
  grad_hess_ret[1] << 12 * xy, 6 * x_sq + 6 * y, 3 * z_sq, 6 * x_sq + 6 * y,
      6 * x, 0, 3 * z_sq, 0, 6 * zx;
  grad_hess_ret[2] << 0, 3 * z_sq, 6 * zy, 3 * z_sq, 0, 6 * zx, 6 * zy, 6 * zx,
      6 * yx;
  return grad_hess_ret;
}

std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >
norm_grad_hess(const Eigen::Matrix<double, Eigen::Dynamic, 1>& inp_vec) {
  using Eigen::Dynamic;
  using Eigen::Matrix;

  std::vector<Matrix<double, Dynamic, Dynamic> > grad_hess;

  for (int i = 0; i < 3; ++i)
    grad_hess.push_back(Matrix<double, Dynamic, Dynamic>(3, 3));
  double y = inp_vec(0);
  double mu = inp_vec(1);
  double sig = inp_vec(2);

  double inv_sigma_cub = 1 / (sig * sig * sig);
  double inv_sigma_four = inv_sigma_cub / sig;
  double y_m_mu = y - mu;
  double norm_113 = 2 * inv_sigma_cub;
  double norm_123 = -norm_113;
  double norm_223 = norm_113;
  double norm_233 = 6 * inv_sigma_four * y_m_mu;
  double norm_133 = -norm_233;
  double norm_333 = norm_123 + 12 * inv_sigma_four / sig * y_m_mu * y_m_mu;

  grad_hess[0] << 0, 0, norm_113, 0, 0, norm_123, norm_113, norm_123, norm_133;

  grad_hess[1] << 0, 0, norm_123, 0, 0, norm_223, norm_123, norm_223, norm_233;

  grad_hess[2] << norm_113, norm_123, norm_133, norm_123, norm_223, norm_233,
      norm_133, norm_233, norm_333;

  return grad_hess;
}

#endif
