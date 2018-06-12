#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
// #include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <stan/math/rev/mat/fun/matrix_exp_multiply.hpp>
#include <stan/math/rev/mat/fun/to_var.hpp>
#include <vector>

template <int N, int M>
inline void test_matrix_exp_multiply_dv() {
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

  // matrix_exp_multiply
  Eigen::Matrix<var, N, M> res_dv = stan::math::matrix_exp_multiply(A, Bv);
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

TEST(MathMatrix, matrix_exp_multiply_dv) {
  test_matrix_exp_multiply_dv<1, 1>();
  test_matrix_exp_multiply_dv<1, 5>();
  test_matrix_exp_multiply_dv<5, 1>();
  test_matrix_exp_multiply_dv<5, 5>();
}

template <int N, int M>
inline void test_matrix_exp_multiply_vd() {
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
  // matrix_exp_multiply
  Eigen::Matrix<var, N, M> res_vd = stan::math::matrix_exp_multiply(Av, B);
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

TEST(MathMatrix, matrix_exp_multiply_vd) {
  test_matrix_exp_multiply_vd<1, 1>();
  test_matrix_exp_multiply_vd<1, 5>();
  test_matrix_exp_multiply_vd<5, 1>();
  test_matrix_exp_multiply_vd<5, 5>();
}

template <int N, int M>
inline void test_matrix_exp_multiply_vv() {
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

  // matrix_exp_multiply
  Eigen::Matrix<var, N, M> res_vv = stan::math::matrix_exp_multiply(Av, Bv);
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

TEST(MathMatrix, matrix_exp_multiply_vv) {
  test_matrix_exp_multiply_vv<1, 1>();
  test_matrix_exp_multiply_vv<1, 5>();
  test_matrix_exp_multiply_vv<5, 1>();
  test_matrix_exp_multiply_vv<5, 5>();
  test_matrix_exp_multiply_vv<10, 2>();
}

// TEST(MathMatrix, matrix_exp_multiply_clock) {
//   using Eigen::MatrixXd;
//   using stan::math::var;
//   using stan::math::value_of;
//   std::srand(2999);

//   for (int i = 1; i < 50; ++i) {
//     MatrixXd A = MatrixXd::Random(i*5, i*5);
//     MatrixXd B = MatrixXd::Random(i*5, 1);

//     auto start = std::chrono::system_clock::now();
//     auto expAB = stan::math::matrix_exp_action(A, B);
//     auto end = std::chrono::system_clock::now();
//     std::chrono::duration<double> elapsed_seconds = end-start;
//     std::cout << i*5 << "\t" << elapsed_seconds.count() << "\t";

//     start = std::chrono::system_clock::now();
//     expAB = stan::math::multiply(stan::math::matrix_exp(A), B);
//     end = std::chrono::system_clock::now();
//     elapsed_seconds = end-start;
//     std::cout << elapsed_seconds.count() << "\n";
//   }
// }

// TEST(MathMatrix, matrix_exp_action_clock_large_norm) {
//   using Eigen::MatrixXd;
//   using stan::math::var;
//   using stan::math::value_of;
//   std::srand(2999);

//   for (int i = 1; i < 50; ++i) {
//     MatrixXd A = 1000 * MatrixXd::Random(i*5, i*5);
//     MatrixXd B = 100 * MatrixXd::Random(i*5, 1);

//     auto start = std::chrono::system_clock::now();
//     auto expAB = stan::math::matrix_exp_action(A, B);
//     auto end = std::chrono::system_clock::now();
//     std::chrono::duration<double> elapsed_seconds = end-start;
//     std::cout << i*5 << "\t" << elapsed_seconds.count() << "\t";

//     start = std::chrono::system_clock::now();
//     expAB = stan::math::multiply(stan::math::matrix_exp(A), B);
//     end = std::chrono::system_clock::now();
//     elapsed_seconds = end-start;
//     std::cout << elapsed_seconds.count() << "\n";
//   }
// }
