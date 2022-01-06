#include <stan/math/mix.hpp>

#include <gtest/gtest.h>
#include <test/unit/util.hpp>

TEST(MathMetaMix, PromoteScalarScalar) {
  using stan::math::promote_scalar_t;
  using stan::math::fvar;
  using stan::math::var;
  using fvar_d = fvar<double>;
  using fvar_v = fvar<var>;
  using complex_d = std::complex<double>;
  using complex_v = std::complex<var>;
  EXPECT_SAME_TYPE(double, promote_scalar_t<double, var>);
  EXPECT_SAME_TYPE(var, promote_scalar_t<var, double>);
  EXPECT_SAME_TYPE(complex_d, promote_scalar_t<complex_d, double>);
  EXPECT_SAME_TYPE(complex_v, promote_scalar_t<complex_v, var>);
}

TEST(MathMetaMix, PromoteScalarMatrix) {
  using stan::math::promote_scalar_t;
  using stan::math::fvar;
  using stan::math::var;
  using fvar_d = fvar<double>;
  using fvar_v = fvar<var>;
  using complex_d = std::complex<double>;
  using complex_v = std::complex<var>;
  using mat_d = Eigen::Matrix<double, -1, -1>;
  using mat_v = Eigen::Matrix<var, -1, -1>;
  using mat_cd = Eigen::Matrix<std::complex<double>, -1, -1>;
  using mat_cv = Eigen::Matrix<std::complex<var>, -1, -1>;
  EXPECT_SAME_TYPE(mat_d, promote_scalar_t<mat_d, double>);
  EXPECT_SAME_TYPE(mat_v, promote_scalar_t<var, mat_d>);
  EXPECT_SAME_TYPE(mat_cd, promote_scalar_t<complex_d, mat_d>);
  EXPECT_SAME_TYPE(mat_cv, promote_scalar_t<complex_v, mat_d>);
}


TEST(MathMetaMix, PromoteScalarStdVector) {
  using stan::math::promote_scalar_t;
  using stan::math::fvar;
  using stan::math::var;
  using fvar_d = fvar<double>;
  using fvar_v = fvar<var>;
  using complex_d = std::complex<double>;
  using complex_v = std::complex<var>;
  using mat_d = Eigen::Matrix<double, -1, -1>;
  using mat_v = Eigen::Matrix<var, -1, -1>;
  using mat_cd = Eigen::Matrix<std::complex<double>, -1, -1>;
  using mat_cv = Eigen::Matrix<std::complex<var>, -1, -1>;
  using std_vec_d = std::vector<double>;
  using std_vec_v = std::vector<var>;
  using std_vec_cd = std::vector<std::complex<double>>;
  using std_vec_cv = std::vector<std::complex<var>>;
  using std_vec_mat_d = std::vector<mat_d>;
  using std_vec_mat_v = std::vector<mat_v>;
  using std_vec_mat_cd = std::vector<mat_cd>;
  using std_vec_mat_cv = std::vector<mat_cv>;

  EXPECT_SAME_TYPE(std_vec_d, promote_scalar_t<std_vec_d, double>);
  EXPECT_SAME_TYPE(std_vec_v, promote_scalar_t<var, std_vec_d>);
  EXPECT_SAME_TYPE(std_vec_cd, promote_scalar_t<complex_d, std_vec_d>);
  EXPECT_SAME_TYPE(std_vec_cv, promote_scalar_t<complex_v, std_vec_d>);
  EXPECT_SAME_TYPE(std_vec_mat_cd, promote_scalar_t<complex_d, std_vec_mat_d>);
  EXPECT_SAME_TYPE(std_vec_mat_cv, promote_scalar_t<complex_v, std_vec_mat_d>);
}
