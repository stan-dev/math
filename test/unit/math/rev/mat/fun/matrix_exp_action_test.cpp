#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
// #include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <stan/math/rev/mat/fun/matrix_exp_action_handler.hpp>
#include <stan/math/rev/mat/fun/matrix_exp_action.hpp>
#include <stan/math/rev/mat/fun/to_var.hpp>


TEST(MathMatrix, matrix_exp_action_diag) {
  using MatrixType = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
  using VectorType = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  stan::math::matrix_exp_action_handler handler;

  {
    double t = 0.5;
    MatrixType m1(2, 2);
    VectorType b(2);
    m1 << 2, 0, 0, 2;
    b << 1, 1;
    auto res = handler.action(m1, b, t);
    EXPECT_NEAR(res(0), M_E, 1.e-8);
    EXPECT_NEAR(res(1), M_E, 1.e-8);
  }

  {
    double t = 1.0;
    MatrixType m1(2, 2);
    VectorType b = VectorType::Random(2);
    m1 << 1, 0, 0, 2;
    auto res = handler.action(m1, b, t);
    EXPECT_NEAR(res(0), b(0)*M_E, 1.e-8);
    EXPECT_NEAR(res(1), b(1)*M_E*M_E, 1.e-8);
  }

  {
    double t = 1.0;
    std::srand(1299);
    MatrixType m1(2, 2);
    VectorType b = VectorType::Random(2);
    m1 << -4.0, 0, 0, -5.0;
    auto res = handler.action(m1, b, t);
    EXPECT_NEAR(res(0), b(0)/(M_E*M_E*M_E*M_E), 1.e-8);
    EXPECT_NEAR(res(1), b(1)/(M_E*M_E*M_E*M_E*M_E), 1.e-8);
  }

  {
    std::srand(999);
    double t = double(std::rand())/RAND_MAX;
    VectorType b = VectorType::Random(5);
    VectorType d = VectorType::Random(5);
    MatrixType m = d.asDiagonal();
    auto res = handler.action(m, b, t);
    for (int i = 0; i < 5; ++i) {
      EXPECT_NEAR(res(i), b(i) * std::exp(t * d(i)), 1.e-8);
    }
  }
}

TEST(MathMatrix, matrix_exp_action_vector) {
  using MatrixType = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
  using VectorType = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  stan::math::matrix_exp_action_handler handler;
  std::srand(999);

  for (size_t n = 2; n < 10; ++n) {
    MatrixType A = MatrixType::Random(n, n);
    VectorType b = VectorType::Random(n);
    VectorType res = handler.action(A, b);
    MatrixType expA = stan::math::matrix_exp(A);
    VectorType expb = expA * b;
    for (size_t i = 0; i < n; ++i) {
      EXPECT_NEAR(res(i), expb(i), 1.e-6);
    }

    int m1, s1, m2, s2;
    const double t1 = 9.9, t2 = 1.0;
    handler.set_approximation_parameter(A, t1, m1, s1);
    A *= t1;
    handler.set_approximation_parameter(A, t2, m2, s2);
    EXPECT_EQ(m1, m2);
    EXPECT_EQ(s1, s2);
  }

}

TEST(MathMatrix, matrix_exp_action_matrix) {
  using MatrixType = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
  stan::math::matrix_exp_action_handler handler;
  std::srand(999);

  constexpr int N = 10;
  constexpr int M = 4;
  Eigen::Matrix<double, N, N> A = Eigen::Matrix<double, N, N>::Random();
  Eigen::Matrix<double, N, M> B = Eigen::Matrix<double, N, M>::Random();

  Eigen::Matrix<double, N, M> res = handler.action(A, B);
  MatrixType Ad(A);
  MatrixType expa = stan::math::matrix_exp(Ad);
  Eigen::Matrix<double, N, M> expb = expa * B;

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < M; ++j) {    
      EXPECT_NEAR(res(i, j), expb(i, j), 1.e-6);          
    }
  }

}

TEST(MathMatrix, matrix_exp_action_matrix_transpose) {
  using MatrixType = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
  stan::math::matrix_exp_action_handler handler;
  std::srand(1999);

  constexpr int N = 10;
  constexpr int M = 4;
  Eigen::Matrix<double, N, N> A = Eigen::Matrix<double, N, N>::Random();
  Eigen::Matrix<double, N, M> B = Eigen::Matrix<double, N, M>::Random();

  Eigen::Matrix<double, N, M> res = handler.action(A.transpose(), B);
  MatrixType Ad(A);
  MatrixType expa = stan::math::matrix_exp(Ad).transpose();
  Eigen::Matrix<double, N, M> expb = expa * B;

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < M; ++j) {    
      EXPECT_NEAR(res(i, j), expb(i, j), 1.e-6);
    }
  }
}


template<int N, int M>
inline
void test_matrix_exp_action_dv() {
  using stan::math::var;
  using stan::math::value_of;
  stan::math::matrix_exp_action_handler handler;
  std::srand(1999);

  Eigen::Matrix<var, N, N> Av = Eigen::Matrix<var, N, N>::Random();
  Eigen::Matrix<var, N, M> Bv = Eigen::Matrix<var, N, M>::Random();
  std::vector<stan::math::var> Bvec(Bv.data(), Bv.data() + Bv.size());
  std::vector<stan::math::var> Avec(Av.data(), Av.data() + Av.size());
  Eigen::Matrix<double, N, N> A = value_of(Av);
  Eigen::Matrix<double, N, M> B = value_of(Bv);

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
  Eigen::Matrix<var, N, M> res_dv = stan::math::matrix_exp_action(A, Bv);

  // compare adjoints
  std::vector<double> g, g0;
  for (int l = 0; l < M; ++l) {
    for (int k = 0; k < N; ++k) {
      stan::math::set_zero_all_adjoints();
      res_dv(k, l).grad(Bvec, g);      
      stan::math::set_zero_all_adjoints();
      expAB(k, l).grad(Bvec, g0);
      for (size_t j = 0; j < g.size(); ++j) {
        EXPECT_NEAR(g[j], g0[j], 1.e-6);
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
    EXPECT_NEAR(g[j], g0[j], 1.e-6);
  }
}

TEST(MathMatrix, matrix_exp_action_dv) {
  test_matrix_exp_action_dv<1, 1>();
  test_matrix_exp_action_dv<1, 5>();
  test_matrix_exp_action_dv<5, 1>();
  test_matrix_exp_action_dv<5, 5>();
}

template<int N, int M>
inline
void test_matrix_exp_action_vd() {
  using stan::math::var;
  using stan::math::value_of;
  stan::math::matrix_exp_action_handler handler;
  std::srand(1999);

  Eigen::Matrix<var, N, N> Av = Eigen::Matrix<var, N, N>::Random();
  Eigen::Matrix<var, N, M> Bv = Eigen::Matrix<var, N, M>::Random();
  std::vector<stan::math::var> Bvec(Bv.data(), Bv.data() + Bv.size());
  std::vector<stan::math::var> Avec(Av.data(), Av.data() + Av.size());
  Eigen::Matrix<double, N, N> A = value_of(Av);
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
  Eigen::Matrix<var, N, M> res_vd = stan::math::matrix_exp_action(Av, B);

  // compare adjoints
  std::vector<double> g, g0;
  for (int l = 0; l < M; ++l) {
    for (int k = 0; k < N; ++k) {
      stan::math::set_zero_all_adjoints();
      res_vd(k, l).grad(Avec, g);      
      stan::math::set_zero_all_adjoints();
      expAB(k, l).grad(Avec, g0);
      for (size_t j = 0; j < g.size(); ++j) {
        EXPECT_NEAR(g[j], g0[j], 1.e-6);
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
    EXPECT_NEAR(g[j], g0[j], 1.e-6);
  }
}

TEST(MathMatrix, matrix_exp_action_vd) {
  test_matrix_exp_action_vd<1, 1>();
  test_matrix_exp_action_vd<1, 5>();
  test_matrix_exp_action_vd<5, 1>();
  test_matrix_exp_action_vd<5, 5>();
}

template<int N, int M>
inline
void test_matrix_exp_action_vv() {
  using stan::math::var;
  using stan::math::value_of;
  stan::math::matrix_exp_action_handler handler;
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
  Eigen::Matrix<var, N, M> res_vv = stan::math::matrix_exp_action(Av, Bv);
  Avec.insert( Avec.end(), Bvec.begin(), Bvec.end());

  // compare adjoints
  std::vector<double> g, g0;
  for (int l = 0; l < M; ++l) {
    for (int k = 0; k < N; ++k) {
      stan::math::set_zero_all_adjoints();
      res_vv(k, l).grad(Avec, g);      
      stan::math::set_zero_all_adjoints();
      expAB(k, l).grad(Avec, g0);
      for (size_t j = 0; j < g.size(); ++j) {
        EXPECT_NEAR(g[j], g0[j], 1.e-6);
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
    EXPECT_NEAR(g[j], g0[j], 1.e-6);
  }
}

TEST(MathMatrix, matrix_exp_action_vv) {
  test_matrix_exp_action_vv<1, 1>();
  test_matrix_exp_action_vv<1, 5>();
  test_matrix_exp_action_vv<5, 1>();
  test_matrix_exp_action_vd<5, 5>();
}
