#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/expect_matrix_eq.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <random>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <vector>

TEST(MathMatrix, matrix_exp_1x1) {
  stan::math::matrix_v m1(1, 1), m2(1, 1), m1_exp;
  Eigen::Matrix<stan::math::var, 1, 1> m3, m3_exp;
  m1 << 0;
  m2 << 1;
  m1_exp = matrix_exp(m1);
  m3 = m1;
  m3_exp = matrix_exp(m3);
  expect_matrix_eq(m2, m1_exp);

  AVEC x = createAVEC(m1(0, 0));
  AVEC y = createAVEC(m3(0, 0));
  VEC g;
  m1_exp(0, 0).grad(x, g);
  EXPECT_FLOAT_EQ(m1_exp(0, 0).val(), g[0]);
  stan::math::set_zero_all_adjoints();
  m3_exp(0, 0).grad(y, g);
  EXPECT_FLOAT_EQ(m3_exp(0, 0).val(), g[0]);
}

TEST(MathMatrix, matrix_exp_2x2) {
  using stan::math::matrix_v;

  // example from Moler & Van Loan, 2003
  for (size_t k = 0; k < 2; k++) {
    for (size_t l = 0; l < 2; l++) {
      AVAR a = -1.0, b = -17.0;

      matrix_v m1(2, 2), m2(2, 2), m1_exp;
      m1 << -2 * a + 3 * b, 1.5 * a - 1.5 * b, -4 * a + 4 * b, 3 * a - 2 * b;
      m2 << -.735759, .551819, -1.471518, 1.103638;
      m1_exp = matrix_exp(m1);
      expect_matrix_eq(m2, m1_exp);

      matrix_v dm1_exp_da(2, 2), dm1_exp_db(2, 2);
      AVAR exp_a = exp(a), exp_b = exp(b);
      dm1_exp_da << -2 * exp_a, 1.5 * exp_a, -4 * exp_a, 3 * exp_a;
      dm1_exp_db << 3 * exp_b, -1.5 * exp_b, 4 * exp_b, -2 * exp_b;

      AVEC x = createAVEC(a, b);
      VEC g;
      m1_exp(k, l).grad(x, g);
      EXPECT_FLOAT_EQ(dm1_exp_da(k, l).val(), g[0]);
      EXPECT_FLOAT_EQ(dm1_exp_db(k, l).val(), g[1]);
    }
  }
}

TEST(MathMatrix, matrix_exp_2x2_2) {
  // make sure matrix_exp doesn't use matrix_exp_2x2,
  // which would return NaN for this matrix
  // Compare to result from http:// comnuan.com/cmnn01015/
  // Don't test derivatives, since goal is to see that
  // matrix_exp picks the right algorithm

  stan::math::matrix_v m(2, 2), exp_m(2, 2);
  m << -0.999984, 0.511211, -0.736924, -0.0826997;

  exp_m << 0.2746483852, 0.2893267425, -0.4170720513, 0.7937977746;

  expect_matrix_eq(exp_m, stan::math::matrix_exp(m));
}

TEST(MathMatrix, matrix_exp_3x3) {
  using stan::math::matrix_v;

  for (size_t k = 0; k < 3; k++) {
    for (size_t l = 0; l < 3; l++) {
      AVAR a = -1.0, b = 2.0, c = 1.0;

      matrix_v m1(3, 3), m2(3, 3), m1_exp;
      m1 << -24 * a + 40 * b - 15 * c, 18 * a - 30 * b + 12 * c,
          5 * a - 8 * b + 3 * c, 20 * b - 20 * c, -15 * b + 16 * c,
          -4 * b + 4 * c, -120 * a + 120 * b, 90 * a - 90 * b, 25 * a - 24 * b;
      m2 << 245.95891, -182.43047, -49.11821, 93.41549, -67.3433, -18.68310,
          842.54120, -631.90590, -168.14036;
      m1_exp = matrix_exp(m1);
      expect_matrix_eq(m2, m1_exp);

      matrix_v dm1_exp_da(3, 3), dm1_exp_db(3, 3), dm1_exp_dc(3, 3);
      AVAR exp_a = exp(a), exp_b = exp(b), exp_c = exp(c);
      dm1_exp_da << -24 * exp_a, 18 * exp_a, 5 * exp_a, 0, 0, 0, -120 * exp_a,
          90 * exp_a, 25 * exp_a;
      dm1_exp_db << 40 * exp_b, -30 * exp_b, -8 * exp_b, 20 * exp_b,
          -15 * exp_b, -4 * exp_b, 120 * exp_b, -90 * exp_b, -24 * exp_b;
      dm1_exp_dc << -15 * exp_c, 12 * exp_c, 3 * exp_c, -20 * exp_c, 16 * exp_c,
          4 * exp_c, 0, 0, 0;

      AVEC x = createAVEC(a, b, c);
      VEC g;
      m1_exp(k, l).grad(x, g);

      double rel_err
          = std::max(dm1_exp_da.cwiseAbs().maxCoeff().val(), std::abs(g[0]))
            * 1e-10;
      EXPECT_NEAR(dm1_exp_da(k, l).val(), g[0], rel_err);

      rel_err = std::max(dm1_exp_db.cwiseAbs().maxCoeff().val(), std::abs(g[1]))
                * 1e-10;
      EXPECT_NEAR(dm1_exp_db(k, l).val(), g[1], rel_err);

      rel_err = std::max(dm1_exp_dc.cwiseAbs().maxCoeff().val(), std::abs(g[2]))
                * 1e-10;
      EXPECT_NEAR(dm1_exp_dc(k, l).val(), g[2], rel_err);
    }
  }
}

TEST(MathMatrix, matrix_exp_25x25) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::matrix_v;
  using stan::math::var;

  int size = 25;

  std::random_device rd;
  std::mt19937 mt(rd());

  // Randomly construct input matrix
  Matrix<double, Dynamic, Dynamic> S = Eigen::MatrixXd::Identity(size, size),
                                   I = Eigen::MatrixXd::Identity(size, size);
  int col1, col2;
  for (int i = 0; i < 5 * size; i++) {
    col1 = rd() % size;
    col2 = rd() % size;
    while (col1 == col2)
      col2 = rd() % size;
    S.col(col1) += S.col(col2) * std::pow(-1, rd());
  }
  Matrix<double, Dynamic, Dynamic> S_inv = stan::math::mdivide_right(I, S);
  Matrix<double, Dynamic, Dynamic> diag_elements(1, size);
  diag_elements.setRandom();

  // Compute matrix exponential (analytically and using function)
  stan::math::start_nested();
  matrix_v diag_elements_var = diag_elements.cast<var>();
  matrix_v A
      = S.cast<var>() * diag_elements_var.asDiagonal() * S_inv.cast<var>();
  matrix_v expm_A = stan::math::matrix_exp(A);
  Matrix<double, Dynamic, Dynamic> exp_diag_elements
      = stan::math::exp(diag_elements);
  Matrix<double, Dynamic, Dynamic> exp_A
      = S * exp_diag_elements.asDiagonal() * S_inv;

  double rel_err = 1e-10
                   * std::max(exp_A.cwiseAbs().maxCoeff(),
                              expm_A.cwiseAbs().maxCoeff().val());

  // Test values
  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++)
      EXPECT_NEAR(exp_A(i, j), expm_A(i, j).val(), rel_err);

  // Test adjoints
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      if (i > 0 || j > 0)
        stan::math::set_zero_all_adjoints_nested();
      grad(expm_A(i, j).vi_);
      std::vector<double> g_abs(size);
      for (int k = 0; k < size; k++)
        g_abs[k] = std::abs(diag_elements_var(k).adj());
      Matrix<double, Dynamic, Dynamic> dA(size, size), dA_exp(size, size);
      for (int k = 0; k < size; k++) {
        dA.setZero();
        dA(k, k) = exp(diag_elements(k));
        dA_exp = S * dA * S_inv;
        rel_err = std::max(dA_exp.cwiseAbs().maxCoeff(), stan::math::max(g_abs))
                  * 1e-10;
        EXPECT_NEAR(dA_exp(i, j), diag_elements_var(k).adj(), rel_err);
      }
    }
  }
}

TEST(MathMatrix, matrix_exp_exceptions) {
  stan::math::matrix_v m1(0, 0), m2(1, 2);
  m2 << 1, 2;

  EXPECT_THROW(matrix_exp(m1), std::invalid_argument);
  EXPECT_THROW(matrix_exp(m2), std::invalid_argument);
}

////////////////////////////////////////////////////////////////
// Overload matrix_exp_action.
////////////////////////////////////////////////////////////////
template <int N, int M>
inline void test_matrix_exp_overload_dv() {
  using stan::math::value_of;
  using stan::math::var;
  std::srand(1999);

  Eigen::Matrix<var, N, N> Av = Eigen::Matrix<var, N, N>::Random();
  Eigen::Matrix<var, N, M> Bv = Eigen::Matrix<var, N, M>::Random();
  std::vector<stan::math::var> Bvec(Bv.data(), Bv.data() + Bv.size());
  std::vector<stan::math::var> Avec(Av.data(), Av.data() + Av.size());
  Eigen::Matrix<double, N, N> A = value_of(Av);

  // brute force
  Eigen::Matrix<double, N, N> expA = stan::math::matrix_exp(A);
  Eigen::Matrix<var, N, M> expAB;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < M; ++j) {
      expAB(i, j) = 0.0;
      for (int k = 0; k < N; ++k) {
        expAB(i, j) += expA(i, k) * Bv(k, j);
      }
    }
  }

  // matrix_exp_action
  Eigen::Matrix<var, N, M> res_dv = stan::math::matrix_exp(A, Bv);
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

TEST(MathMatrix, matrix_exp_overload_dv) {
  test_matrix_exp_overload_dv<1, 1>();
  test_matrix_exp_overload_dv<1, 5>();
  test_matrix_exp_overload_dv<5, 1>();
  test_matrix_exp_overload_dv<5, 5>();
}

template <int N, int M>
inline void test_matrix_exp_overload_vd() {
  using stan::math::value_of;
  using stan::math::var;
  std::srand(1999);

  Eigen::Matrix<var, N, N> Av = Eigen::Matrix<var, N, N>::Random();
  Eigen::Matrix<var, N, M> Bv = Eigen::Matrix<var, N, M>::Random();
  std::vector<stan::math::var> Bvec(Bv.data(), Bv.data() + Bv.size());
  std::vector<stan::math::var> Avec(Av.data(), Av.data() + Av.size());
  Eigen::Matrix<double, N, M> B = value_of(Bv);

  // brute force
  Eigen::Matrix<var, N, N> expA = stan::math::matrix_exp(Av);
  Eigen::Matrix<var, N, M> expAB;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < M; ++j) {
      expAB(i, j) = 0.0;
      for (int k = 0; k < N; ++k) {
        expAB(i, j) += expA(i, k) * B(k, j);
      }
    }
  }
  // matrix_exp_action
  Eigen::Matrix<var, N, M> res_vd = stan::math::matrix_exp(Av, B);
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

TEST(MathMatrix, matrix_exp_overload_vd) {
  test_matrix_exp_overload_vd<1, 1>();
  test_matrix_exp_overload_vd<1, 5>();
  test_matrix_exp_overload_vd<5, 1>();
  test_matrix_exp_overload_vd<5, 5>();
}

template <int N, int M>
inline void test_matrix_exp_overload_vv() {
  using stan::math::value_of;
  using stan::math::var;
  std::srand(2999);

  Eigen::Matrix<var, N, N> Av = Eigen::Matrix<var, N, N>::Random();
  Eigen::Matrix<var, N, M> Bv = Eigen::Matrix<var, N, M>::Random();
  std::vector<stan::math::var> Bvec(Bv.data(), Bv.data() + Bv.size());
  std::vector<stan::math::var> Avec(Av.data(), Av.data() + Av.size());

  // brute force
  Eigen::Matrix<var, N, N> expA = stan::math::matrix_exp(Av);
  Eigen::Matrix<var, N, M> expAB;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < M; ++j) {
      expAB(i, j) = 0.0;
      for (int k = 0; k < N; ++k) {
        expAB(i, j) += expA(i, k) * Bv(k, j);
      }
    }
  }

  // matrix_exp_action
  Eigen::Matrix<var, N, M> res_vv = stan::math::matrix_exp(Av, Bv);
  for (int l = 0; l < res_vv.size(); ++l) {
    EXPECT_FLOAT_EQ(res_vv(l).val(), expAB(l).val());
  }
  Avec.insert(Avec.end(), Bvec.begin(), Bvec.end());
  // compare adjoints
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

TEST(MathMatrix, matrix_exp_overload_vv) {
  test_matrix_exp_overload_vv<1, 1>();
  test_matrix_exp_overload_vv<1, 5>();
  test_matrix_exp_overload_vv<5, 1>();
  test_matrix_exp_overload_vv<5, 5>();
  test_matrix_exp_overload_vv<10, 2>();
}
