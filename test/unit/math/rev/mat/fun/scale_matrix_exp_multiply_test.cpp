#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
// #include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <stan/math/rev/mat/fun/scale_matrix_exp_multiply.hpp>
#include <stan/math/rev/mat/fun/to_var.hpp>
#include <vector>

inline void test_scale_matrix_exp_multiply_dv(int N, int M) {
  using stan::math::value_of;
  using stan::math::var;

  std::srand(1999);

  Eigen::Matrix<var, -1, -1> Av = Eigen::Matrix<var, -1, -1>::Random(N, N);
  Eigen::Matrix<var, -1, -1> Bv = Eigen::Matrix<var, -1, -1>::Random(N, M);
  std::vector<stan::math::var> Avec = stan::math::to_array_1d(Av);
  std::vector<stan::math::var> Bvec = stan::math::to_array_1d(Bv);
  Eigen::MatrixXd A = value_of(Av);

  // brute force
  Eigen::Matrix<var, -1, -1> expAB
      = stan::math::multiply(stan::math::matrix_exp(A), Bv);

  // matrix_exp_multiply
  const double t = 1.0;
  Eigen::Matrix<var, -1, -1> res_dv
      = stan::math::scale_matrix_exp_multiply(t, A, Bv);
  EXPECT_EQ(res_dv.size(), expAB.size());
  for (int l = 0; l < res_dv.size(); ++l) {
    EXPECT_FLOAT_EQ(res_dv(l).val(), expAB(l).val());
  }

  // compare adjoints
  std::vector<double> g, g0;
  for (int l = 0; l < M; ++l) {
    for (int k = 0; k < N; ++k) {
      stan::math::set_zero_all_adjoints();
      res_dv(k, l).grad(Bvec, g);
      stan::math::set_zero_all_adjoints();
      expAB(k, l).grad(Bvec, g0);
      for (size_t j = 0; j < g.size(); ++j) {
        EXPECT_FLOAT_EQ(g[j], g0[j]);
      }
    }
  }

  // test a single function of expA*B
  var f = sum(res_dv);
  var f0 = sum(expAB);
  stan::math::set_zero_all_adjoints();
  f.grad(Bvec, g);
  stan::math::set_zero_all_adjoints();
  f0.grad(Bvec, g0);
  for (size_t j = 0; j < g.size(); ++j) {
    EXPECT_FLOAT_EQ(g[j], g0[j]);
  }
}

TEST(MathMatrix, scale_matrix_exp_multiply_dv) {
  test_scale_matrix_exp_multiply_dv(1, 1);
  test_scale_matrix_exp_multiply_dv(1, 5);
  test_scale_matrix_exp_multiply_dv(5, 1);
  test_scale_matrix_exp_multiply_dv(5, 5);
}

inline void test_scale_matrix_exp_multiply_vd(int N, int M) {
  using stan::math::value_of;
  using stan::math::var;

  std::srand(1999);

  Eigen::Matrix<var, -1, -1> Av = Eigen::Matrix<var, -1, -1>::Random(N, N);
  Eigen::Matrix<var, -1, -1> Bv = Eigen::Matrix<var, -1, -1>::Random(N, M);
  std::vector<stan::math::var> Avec = stan::math::to_array_1d(Av);
  std::vector<stan::math::var> Bvec = stan::math::to_array_1d(Bv);
  Eigen::MatrixXd B = value_of(Bv);

  // brute force
  Eigen::Matrix<var, -1, -1> expAB
      = stan::math::multiply(stan::math::matrix_exp(Av), B);

  // matrix_exp_multiply
  const double t = 1.0;
  Eigen::Matrix<var, -1, -1> res_vd
      = stan::math::scale_matrix_exp_multiply(t, Av, B);
  EXPECT_EQ(res_vd.size(), expAB.size());
  for (int l = 0; l < res_vd.size(); ++l) {
    EXPECT_FLOAT_EQ(res_vd(l).val(), expAB(l).val());
  }

  // compare adjoints
  std::vector<double> g, g0;
  for (int l = 0; l < M; ++l) {
    for (int k = 0; k < N; ++k) {
      stan::math::set_zero_all_adjoints();
      res_vd(k, l).grad(Avec, g);
      stan::math::set_zero_all_adjoints();
      expAB(k, l).grad(Avec, g0);
      for (size_t j = 0; j < g.size(); ++j) {
        EXPECT_FLOAT_EQ(g[j], g0[j]);
      }
    }
  }

  // test a single function of expA*B
  var f = sum(res_vd);
  var f0 = sum(expAB);
  stan::math::set_zero_all_adjoints();
  f.grad(Avec, g);
  stan::math::set_zero_all_adjoints();
  f0.grad(Avec, g0);
  for (size_t j = 0; j < g.size(); ++j) {
    EXPECT_FLOAT_EQ(g[j], g0[j]);
  }
}

TEST(MathMatrix, scale_matrix_exp_multiply_vd) {
  test_scale_matrix_exp_multiply_vd(1, 1);
  test_scale_matrix_exp_multiply_vd(1, 5);
  test_scale_matrix_exp_multiply_vd(5, 1);
  test_scale_matrix_exp_multiply_vd(5, 5);
}

inline void test_scale_matrix_exp_multiply_vv(int N, int M) {
  using stan::math::value_of;
  using stan::math::var;

  std::srand(2999);

  Eigen::Matrix<var, -1, -1> Av = Eigen::Matrix<var, -1, -1>::Random(N, N);
  Eigen::Matrix<var, -1, -1> Bv = Eigen::Matrix<var, -1, -1>::Random(N, M);
  std::vector<stan::math::var> Avec = stan::math::to_array_1d(Av);
  std::vector<stan::math::var> Bvec = stan::math::to_array_1d(Bv);

  // brute force
  Eigen::Matrix<var, -1, -1> expAB
      = stan::math::multiply(stan::math::matrix_exp(Av), Bv);

  // matrix_exp_multiply
  const double t = 1.0;
  Eigen::Matrix<var, -1, -1> res_vv
      = stan::math::scale_matrix_exp_multiply(t, Av, Bv);
  EXPECT_EQ(res_vv.size(), expAB.size());
  for (int l = 0; l < res_vv.size(); ++l) {
    EXPECT_FLOAT_EQ(res_vv(l).val(), expAB(l).val());
  }

  // compare adjoints
  Avec.insert(Avec.end(), Bvec.begin(), Bvec.end());
  std::vector<double> g, g0;
  for (int l = 0; l < M; ++l) {
    for (int k = 0; k < N; ++k) {
      stan::math::set_zero_all_adjoints();
      res_vv(k, l).grad(Avec, g);
      stan::math::set_zero_all_adjoints();
      expAB(k, l).grad(Avec, g0);
      for (size_t j = 0; j < g.size(); ++j) {
        EXPECT_FLOAT_EQ(g[j], g0[j]);
      }
    }
  }

  // test a single function of expA*B
  var f = sum(res_vv);
  var f0 = sum(expAB);
  stan::math::set_zero_all_adjoints();
  f.grad(Avec, g);
  stan::math::set_zero_all_adjoints();
  f0.grad(Avec, g0);
  for (size_t j = 0; j < g.size(); ++j) {
    EXPECT_FLOAT_EQ(g[j], g0[j]);
  }
}

TEST(MathMatrix, scale_matrix_exp_multiply_vv) {
  test_scale_matrix_exp_multiply_vv(1, 1);
  test_scale_matrix_exp_multiply_vv(1, 5);
  test_scale_matrix_exp_multiply_vv(5, 1);
  test_scale_matrix_exp_multiply_vv(5, 5);
  test_scale_matrix_exp_multiply_vv(8, 2);
}

TEST(MathMatrix, scale_matrix_exp_multiply_exception) {
  using stan::math::var;
  const double t = 1.0;
  {  // nonzero size
    Eigen::Matrix<var, -1, -1> A(0, 0);
    Eigen::Matrix<var, -1, -1> B = Eigen::Matrix<var, -1, -1>::Random(1, 2);
    EXPECT_THROW(scale_matrix_exp_multiply(t, A, B), std::invalid_argument);
    EXPECT_THROW(scale_matrix_exp_multiply(t, B, A), std::invalid_argument);
  }

  {  // multiplicable
    Eigen::Matrix<var, -1, -1> A = Eigen::Matrix<var, -1, -1>::Random(2, 2);
    Eigen::Matrix<var, -1, -1> B = Eigen::Matrix<var, -1, -1>::Random(3, 2);
    EXPECT_THROW(scale_matrix_exp_multiply(t, A, B), std::invalid_argument);
  }

  {  // square
    Eigen::Matrix<var, -1, -1> A = Eigen::Matrix<var, -1, -1>::Random(2, 3);
    Eigen::Matrix<var, -1, -1> B = Eigen::Matrix<var, -1, -1>::Random(3, 2);
    EXPECT_THROW(scale_matrix_exp_multiply(t, A, B), std::invalid_argument);
  }
}
