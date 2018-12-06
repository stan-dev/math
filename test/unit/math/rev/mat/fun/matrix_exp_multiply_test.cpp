#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
// #include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <stan/math/rev/mat/fun/matrix_exp_multiply.hpp>
#include <stan/math/rev/mat/fun/to_var.hpp>
#include <vector>

TEST(MathRev, matrix_multiply_exp__vd__segfault) {
  Eigen::Matrix<stan::math::var, -1, -1> A(5, 5);
  A << -0.96871, 0.398827, 0.241306, 0.741373, 0.108926, 0.888077, -0.915624,
      -0.373344, 0.255238, 0.717304, -0.0899219, -0.898862, -0.800546,
      -0.222652, -0.271382, 0.683227, 0.827031, -0.780702, -0.104228, 0.885106,
      -0.996585, -0.097802, 0.739617, 0.235266, -0.0247717;

  Eigen::Matrix<double, -1, -1> B(5, 5);
  B << -0.96871, 0.398827, 0.241306, 0.741373, 0.108926, 0.888077, -0.915624,
      -0.373344, 0.255238, 0.717304, -0.0899219, -0.898862, -0.800546,
      -0.222652, -0.271382, 0.683227, 0.827031, -0.780702, -0.104228, 0.885106,
      -0.996585, -0.097802, 0.739617, 0.235266, -0.0247717;

  Eigen::Matrix<stan::math::var, -1, -1> result
      = stan::math::matrix_exp_multiply(A, B);
  // use this to debug how large the stack is:
  // stan::math::print_stack(std::cout);

  std::vector<double> gradients;

  std::vector<stan::math::var> vars = stan::math::to_array_1d(A);
  result(0, 0).grad(vars, gradients);

  ASSERT_EQ(25, gradients.size());
  EXPECT_FLOAT_EQ(gradients[0], -0.318076261);
  EXPECT_FLOAT_EQ(gradients[1], -0.103184184);
  EXPECT_FLOAT_EQ(gradients[2], -0.00582563332);
  EXPECT_FLOAT_EQ(gradients[3], -0.178361665);
  EXPECT_FLOAT_EQ(gradients[4], -0.111540085);
  EXPECT_FLOAT_EQ(gradients[5], 0.131275181);
  EXPECT_FLOAT_EQ(gradients[6], 0.068209884);
  EXPECT_FLOAT_EQ(gradients[7], 0.000694387545);
  EXPECT_FLOAT_EQ(gradients[8], 0.116771912);
  EXPECT_FLOAT_EQ(gradients[9], 0.0798794545);
  EXPECT_FLOAT_EQ(gradients[10], -0.0871212464);
  EXPECT_FLOAT_EQ(gradients[11], -0.0246279672);
  EXPECT_FLOAT_EQ(gradients[12], -0.00226500129);
  EXPECT_FLOAT_EQ(gradients[13], -0.0429290518);
  EXPECT_FLOAT_EQ(gradients[14], -0.0246780963);
  EXPECT_FLOAT_EQ(gradients[15], 0.205020349);
  EXPECT_FLOAT_EQ(gradients[16], 0.075942319);
  EXPECT_FLOAT_EQ(gradients[17], 0.00365725736);
  EXPECT_FLOAT_EQ(gradients[18], 0.131097209);
  EXPECT_FLOAT_EQ(gradients[19], 0.0830377219);
  EXPECT_FLOAT_EQ(gradients[20], -0.477182156);
  EXPECT_FLOAT_EQ(gradients[21], -0.133762948);
  EXPECT_FLOAT_EQ(gradients[22], -0.0104704174);
  EXPECT_FLOAT_EQ(gradients[23], -0.232294833);
  EXPECT_FLOAT_EQ(gradients[24], -0.13877342);

  stan::math::recover_memory();
}

inline void test_matrix_exp_multiply_dv(int N, int M) {
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
  Eigen::Matrix<var, -1, -1> res_dv = stan::math::matrix_exp_multiply(A, Bv);
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

TEST(MathMatrix, matrix_exp_multiply_dv) {
  test_matrix_exp_multiply_dv(1, 1);
  test_matrix_exp_multiply_dv(1, 5);
  test_matrix_exp_multiply_dv(5, 1);
  test_matrix_exp_multiply_dv(5, 5);
}

inline void test_matrix_exp_multiply_vd(int N, int M) {
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
  Eigen::Matrix<var, -1, -1> res_vd = stan::math::matrix_exp_multiply(Av, B);
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

TEST(MathMatrix, matrix_exp_multiply_vd) {
  test_matrix_exp_multiply_vd(1, 1);
  test_matrix_exp_multiply_vd(1, 5);
  test_matrix_exp_multiply_vd(5, 1);
  test_matrix_exp_multiply_vd(5, 5);
}

inline void test_matrix_exp_multiply_vv(int N, int M) {
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
  Eigen::Matrix<var, -1, -1> res_vv = stan::math::matrix_exp_multiply(Av, Bv);
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

TEST(MathMatrix, matrix_exp_multiply_vv) {
  test_matrix_exp_multiply_vv(1, 1);
  test_matrix_exp_multiply_vv(1, 5);
  test_matrix_exp_multiply_vv(5, 1);
  test_matrix_exp_multiply_vv(5, 5);
  test_matrix_exp_multiply_vv(8, 2);
}

TEST(MathMatrix, matrix_exp_multiply_exception) {
  using stan::math::matrix_exp_multiply;
  using stan::math::var;
  {  // nonzero size
    Eigen::Matrix<var, -1, -1> A(0, 0);
    Eigen::Matrix<var, -1, -1> B = Eigen::Matrix<var, -1, -1>::Random(1, 2);
    EXPECT_THROW(matrix_exp_multiply(A, B), std::invalid_argument);
    EXPECT_THROW(matrix_exp_multiply(B, A), std::invalid_argument);
  }

  {  // multiplicable
    Eigen::Matrix<var, -1, -1> A = Eigen::Matrix<var, -1, -1>::Random(2, 2);
    Eigen::Matrix<var, -1, -1> B = Eigen::Matrix<var, -1, -1>::Random(3, 2);
    EXPECT_THROW(matrix_exp_multiply(A, B), std::invalid_argument);
  }

  {  // square
    Eigen::Matrix<var, -1, -1> A = Eigen::Matrix<var, -1, -1>::Random(2, 3);
    Eigen::Matrix<var, -1, -1> B = Eigen::Matrix<var, -1, -1>::Random(3, 2);
    EXPECT_THROW(matrix_exp_multiply(A, B), std::invalid_argument);
  }
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
