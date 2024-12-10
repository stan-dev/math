#ifndef STAN_MATH_FWD_FUN_LOG_ADD_EXP_HPP
#define STAN_MATH_FWD_FUN_LOG_ADD_EXP_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/log_add_exp.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

// Overload for fvar and fvar
template <typename T>
inline fvar<T> log_add_exp(const fvar<T>& x1, const fvar<T>& x2) {
  auto val = stan::math::log_add_exp(x1.val_, x2.val_);

  auto exp_x1 = stan::math::exp(x1.val_);
  auto exp_x2 = stan::math::exp(x2.val_);
  auto sum_exp = exp_x1 + exp_x2;

  auto grad1 = exp_x1 / sum_exp;
  auto grad2 = exp_x2 / sum_exp;

  return fvar<T>(val, x1.d_ * grad1 + x2.d_ * grad2);
}

template <typename T>
inline fvar<T> log_add_exp(const fvar<T>& x1, double x2) {
  if (x1.val_ == NEGATIVE_INFTY) {
    return fvar<T>(x2, 0.0);  // log_add_exp(-∞, b) = b
  }
  return log_add_exp(x2, x1);
}

template <typename T>
inline fvar<T> log_add_exp(double x1, const fvar<T>& x2) {
  if (x2.val_ == NEGATIVE_INFTY) {
    return fvar<T>(x1, 0.0);  // log_add_exp(a, -∞) = a
  }
  auto val = stan::math::log_add_exp(x1, x2.val_);
  auto exp_x2 = stan::math::exp(x2.val_);
  auto grad = exp_x2 / (stan::math::exp(x1) + exp_x2);
  return fvar<T>(val, x2.d_ * grad);
}

// Overload for matrices of fvar
template <typename T>
inline Eigen::Matrix<fvar<T>, -1, -1> log_add_exp(
    const Eigen::Matrix<fvar<T>, -1, -1>& a,
    const Eigen::Matrix<fvar<T>, -1, -1>& b) {
  using fvar_mat_type = Eigen::Matrix<fvar<T>, -1, -1>;
  fvar_mat_type result(a.rows(), a.cols());

  // Check for empty inputs
  if (a.size() == 0 || b.size() == 0) {
    throw std::invalid_argument("Input containers must not be empty.");
  }

  // Check for NaN
  if (a.array().isNaN().any() || b.array().isNaN().any()) {
    result.setConstant(fvar<T>(std::numeric_limits<double>::quiet_NaN()));
    return result;
  }

  // Check for infinity
  if (a.array().isInf().any() || b.array().isInf().any()) {
    result.setConstant(fvar<T>(std::numeric_limits<double>::quiet_NaN()));
    return result;
  }

  // Apply the log_add_exp operation directly
  for (int i = 0; i < a.rows(); ++i) {
    for (int j = 0; j < a.cols(); ++j) {
      result(i, j) = stan::math::log_add_exp(a(i, j), b(i, j));
    }
  }

  return result;  // Return the result matrix
}

// Overload for Eigen vectors
template <typename T>
inline Eigen::Matrix<fvar<T>, -1, 1> log_add_exp(
    const Eigen::Matrix<fvar<T>, -1, 1>& a,
    const Eigen::Matrix<fvar<T>, -1, 1>& b) {
  using fvar_vec_type = Eigen::Matrix<fvar<T>, -1, 1>;
  fvar_vec_type result(a.rows());

  // Check for empty inputs
  if (a.size() == 0 || b.size() == 0) {
    throw std::invalid_argument("Input containers must not be empty.");
  }

  // Check for NaN
  if (a.array().isNaN().any() || b.array().isNaN().any()) {
    result.setConstant(fvar<T>(std::numeric_limits<double>::quiet_NaN()));
    return result;
  }

  // Check for infinity
  if (a.array().isInf().any() || b.array().isInf().any()) {
    result.setConstant(fvar<T>(std::numeric_limits<double>::quiet_NaN()));
    return result;
  }

  // Apply the log_add_exp operation directly
  for (int i = 0; i < a.rows(); ++i) {
    result(i) = stan::math::log_add_exp(a(i), b(i));
  }

  return result;  // Return the result vector
}

// Specialization for nested fvar types
template <typename T>
inline auto log_add_exp(
    const Eigen::Matrix<stan::math::fvar<stan::math::fvar<double>>, -1, -1>& a,
    const Eigen::Matrix<stan::math::fvar<stan::math::fvar<double>>, -1, -1>&
        b) {
  using nested_fvar_mat_type
      = Eigen::Matrix<stan::math::fvar<stan::math::fvar<double>>, -1, -1>;
  nested_fvar_mat_type result(a.rows(), a.cols());

  // Check for empty inputs
  if (a.size() == 0 || b.size() == 0) {
    throw std::invalid_argument("Input containers must not be empty.");
  }

  // Check for NaN
  if (a.array().isNaN().any() || b.array().isNaN().any()) {
    result.setConstant(stan::math::fvar<stan::math::fvar<double>>(
        std::numeric_limits<double>::quiet_NaN()));
    return result;
  }

  // Check for infinity
  if (a.array().isInf().any() || b.array().isInf().any()) {
    result.setConstant(stan::math::fvar<stan::math::fvar<double>>(
        std::numeric_limits<double>::quiet_NaN()));
    return result;
  }

  // Implement the logic for log_add_exp for nested fvar types
  for (int i = 0; i < a.rows(); ++i) {
    for (int j = 0; j < a.cols(); ++j) {
      auto inner_a = a(i, j);
      auto inner_b = b(i, j);
      result(i, j) = stan::math::log_add_exp(inner_a, inner_b);
    }
  }

  return result;  // Return the result matrix
}

}  // namespace math
}  // namespace stan

#endif
