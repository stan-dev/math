#ifndef TEST_UNIT_MATH_MIX_FUN_MULTIPLY_UTIL_HPP
#define TEST_UNIT_MATH_MIX_FUN_MULTIPLY_UTIL_HPP
#include <test/unit/math/test_ad.hpp>

template <typename T>
void instantiate_multiply() {
  using stan::math::multiply;
  Eigen::Matrix<double, -1, -1> d_mat(2, 2);
  d_mat << 1, 2, 3, 4;
  Eigen::Matrix<T, -1, -1> v_mat(2, 2);
  v_mat << 1, 2, 3, 4;
  Eigen::Matrix<std::complex<double>, -1, -1> cd_mat(2, 2);
  cd_mat << 1, 2, 3, 4;
  Eigen::Matrix<std::complex<T>, -1, -1> cv_mat(2, 2);
  cv_mat << 1, 2, 3, 4;

  auto d_d_mat = stan::math::eval(multiply(d_mat, d_mat));
  auto d_v_mat = stan::math::eval(multiply(d_mat, v_mat));
  auto d_cd_mat = stan::math::eval(multiply(d_mat, cd_mat));
  auto d_cv_mat = stan::math::eval(multiply(d_mat, cv_mat));

  auto v_d_mat = stan::math::eval(multiply(v_mat, d_mat));
  auto v_v_mat = stan::math::eval(multiply(v_mat, v_mat));
  auto v_cd_mat = stan::math::eval(multiply(v_mat, cd_mat));
  auto v_cv_mat = stan::math::eval(multiply(v_mat, cv_mat));

  auto cd_d_mat = stan::math::eval(multiply(cd_mat, d_mat));
  auto cd_v_mat = stan::math::eval(multiply(cd_mat, v_mat));
  auto cd_cd_mat = stan::math::eval(multiply(cd_mat, cd_mat));
  auto cd_cv_mat = stan::math::eval(multiply(cd_mat, cv_mat));

  auto cv_d_mat = stan::math::eval(multiply(cv_mat, d_mat));
  auto cv_v_mat = stan::math::eval(multiply(cv_mat, v_mat));
  auto cv_cd_mat = stan::math::eval(multiply(cv_mat, cd_mat));
  auto cv_cv_mat = stan::math::eval(multiply(cv_mat, cv_mat));

  Eigen::Matrix<double, -1, 1> d_vec(2);
  d_vec << 1, 2;
  Eigen::Matrix<T, -1, 1> v_vec(2);
  v_vec << 1, 2;
  Eigen::Matrix<std::complex<double>, -1, 1> cd_vec(2);
  cd_vec << 1, 2;
  Eigen::Matrix<std::complex<T>, -1, 1> cv_vec(2);
  cv_vec << 1, 2;

  auto d_d_vec_mat = stan::math::eval(multiply(d_mat, d_vec));
  auto d_v_vec_mat = stan::math::eval(multiply(d_mat, v_vec));
  auto d_cd_vec_mat = stan::math::eval(multiply(d_mat, cd_vec));
  auto d_cv_vec_mat = stan::math::eval(multiply(d_mat, cv_vec));

  auto v_d_vec_mat = stan::math::eval(multiply(v_mat, d_vec));
  auto v_v_vec_mat = stan::math::eval(multiply(v_mat, v_vec));
  auto v_cd_vec_mat = stan::math::eval(multiply(v_mat, cd_vec));
  auto v_cv_vec_mat = stan::math::eval(multiply(v_mat, cv_vec));

  auto cd_d_vec_mat = stan::math::eval(multiply(cd_mat, d_vec));
  auto cd_v_vec_mat = stan::math::eval(multiply(cd_mat, v_vec));
  auto cd_cd_vec_mat = stan::math::eval(multiply(cd_mat, cd_vec));
  auto cd_cv_vec_mat = stan::math::eval(multiply(cd_mat, cv_vec));

  auto cv_d_vec_mat = stan::math::eval(multiply(cv_mat, d_vec));
  auto cv_v_vec_mat = stan::math::eval(multiply(cv_mat, v_vec));
  auto cv_cd_vec_mat = stan::math::eval(multiply(cv_mat, cd_vec));
  auto cv_cv_vec_mat = stan::math::eval(multiply(cv_mat, cv_vec));

  Eigen::Matrix<double, 1, -1> d_rowvec(2);
  d_rowvec << 1, 2;
  Eigen::Matrix<T, 1, -1> v_rowvec(2);
  v_rowvec << 1, 2;
  Eigen::Matrix<std::complex<double>, 1, -1> cd_rowvec(2);
  cd_rowvec << 1, 2;
  Eigen::Matrix<std::complex<T>, 1, -1> cv_rowvec(2);
  cv_rowvec << 1, 2;

  auto d_d_dot_prod = stan::math::eval(multiply(d_vec, d_rowvec));
  auto d_v_dot_prod = stan::math::eval(multiply(d_vec, v_rowvec));
  auto d_cd_dot_prod = stan::math::eval(multiply(d_vec, cd_rowvec));
  auto d_cv_dot_prod = stan::math::eval(multiply(d_vec, cv_rowvec));

  auto v_d_dot_prod = stan::math::eval(multiply(v_vec, d_rowvec));
  auto v_v_dot_prod = stan::math::eval(multiply(v_vec, v_rowvec));
  auto v_cd_dot_prod = stan::math::eval(multiply(v_vec, cd_rowvec));
  auto v_cv_dot_prod = stan::math::eval(multiply(v_vec, cv_rowvec));

  auto cd_d_dot_prod = stan::math::eval(multiply(cd_vec, d_rowvec));
  auto cd_v_dot_prod = stan::math::eval(multiply(cd_vec, v_rowvec));
  auto cd_cd_dot_prod = stan::math::eval(multiply(cd_vec, cd_rowvec));
  auto cd_cv_dot_prod = stan::math::eval(multiply(cd_vec, cv_rowvec));

  auto cv_d_dot_prod = stan::math::eval(multiply(cv_vec, d_rowvec));
  auto cv_v_dot_prod = stan::math::eval(multiply(cv_vec, v_rowvec));
  auto cv_cd_dot_prod = stan::math::eval(multiply(cv_vec, cd_rowvec));
  auto cv_cv_dot_prod = stan::math::eval(multiply(cv_vec, cv_rowvec));

  auto d_d_outer_prod = stan::math::eval(multiply(d_rowvec, d_vec));
  auto d_v_outer_prod = stan::math::eval(multiply(d_rowvec, v_vec));
  auto d_cd_outer_prod = stan::math::eval(multiply(d_rowvec, cd_vec));
  auto d_cv_outer_prod = stan::math::eval(multiply(d_rowvec, cv_vec));

  auto v_d_outer_prod = stan::math::eval(multiply(v_rowvec, d_vec));
  auto v_v_outer_prod = stan::math::eval(multiply(v_rowvec, v_vec));
  auto v_cd_outer_prod = stan::math::eval(multiply(v_rowvec, cd_vec));
  auto v_cv_outer_prod = stan::math::eval(multiply(v_rowvec, cv_vec));

  auto cd_d_outer_prod = stan::math::eval(multiply(cd_rowvec, d_vec));
  auto cd_v_outer_prod = stan::math::eval(multiply(cd_rowvec, v_vec));
  auto cd_cd_outer_prod = stan::math::eval(multiply(cd_rowvec, cd_vec));
  auto cd_cv_outer_prod = stan::math::eval(multiply(cd_rowvec, cv_vec));

  auto cv_d_outer_prod = stan::math::eval(multiply(cv_rowvec, d_vec));
  auto cv_v_outer_prod = stan::math::eval(multiply(cv_rowvec, v_vec));
  auto cv_cd_outer_prod = stan::math::eval(multiply(cv_rowvec, cd_vec));
  auto cv_cv_outer_prod = stan::math::eval(multiply(cv_rowvec, cv_vec));

  auto d_d_rowvec_mat = stan::math::eval(multiply(d_rowvec, d_mat));
  auto d_v_rowvec_mat = stan::math::eval(multiply(d_rowvec, v_mat));
  auto d_cd_rowvec_mat = stan::math::eval(multiply(d_rowvec, cd_mat));
  auto d_cv_rowvec_mat = stan::math::eval(multiply(d_rowvec, cv_mat));

  auto v_d_rowvec_mat = stan::math::eval(multiply(v_rowvec, d_mat));
  auto v_v_rowvec_mat = stan::math::eval(multiply(v_rowvec, v_mat));
  auto v_cd_rowvec_mat = stan::math::eval(multiply(v_rowvec, cd_mat));
  auto v_cv_rowvec_mat = stan::math::eval(multiply(v_rowvec, cv_mat));

  auto cd_d_rowvec_mat = stan::math::eval(multiply(cd_rowvec, d_mat));
  auto cd_v_rowvec_mat = stan::math::eval(multiply(cd_rowvec, v_mat));
  auto cd_cd_rowvec_mat = stan::math::eval(multiply(cd_rowvec, cd_mat));
  auto cd_cv_rowvec_mat = stan::math::eval(multiply(cd_rowvec, cv_mat));

  auto cv_d_rowvec_mat = stan::math::eval(multiply(cv_rowvec, d_mat));
  auto cv_v_rowvec_mat = stan::math::eval(multiply(cv_rowvec, v_mat));
  auto cv_cd_rowvec_mat = stan::math::eval(multiply(cv_rowvec, cd_mat));
  auto cv_cv_rowvec_mat = stan::math::eval(multiply(cv_rowvec, cv_mat));

  auto d_d_mat_vec = stan::math::eval(multiply(d_mat, d_vec));
  auto d_v_mat_vec = stan::math::eval(multiply(d_mat, v_vec));
  auto d_cd_mat_vec = stan::math::eval(multiply(d_mat, cd_vec));
  auto d_cv_mat_vec = stan::math::eval(multiply(d_mat, cv_vec));

  auto v_d_mat_vec = stan::math::eval(multiply(v_mat, d_vec));
  auto v_v_mat_vec = stan::math::eval(multiply(v_mat, v_vec));
  auto v_cd_mat_vec = stan::math::eval(multiply(v_mat, cd_vec));
  auto v_cv_mat_vec = stan::math::eval(multiply(v_mat, cv_vec));

  auto cd_d_mat_vec = stan::math::eval(multiply(cd_mat, d_vec));
  auto cd_v_mat_vec = stan::math::eval(multiply(cd_mat, v_vec));
  auto cd_cd_mat_vec = stan::math::eval(multiply(cd_mat, cd_vec));
  auto cd_cv_mat_vec = stan::math::eval(multiply(cd_mat, cv_vec));

  auto cv_d_mat_vec = stan::math::eval(multiply(cv_mat, d_vec));
  auto cv_v_mat_vec = stan::math::eval(multiply(cv_mat, v_vec));
  auto cv_cd_mat_vec = stan::math::eval(multiply(cv_mat, cd_vec));
  auto cv_cv_mat_vec = stan::math::eval(multiply(cv_mat, cv_vec));
}

#endif
