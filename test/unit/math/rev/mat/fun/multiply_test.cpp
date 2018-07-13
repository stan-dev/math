#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

// multiply operates on two matrices A X B
// A (n, m)
// B (m, k)
// If stacked col-wise like vec operator
// first n*m elements are for A
// second m*k elements starting with n*m + 1-th element (indexed by
// n*m (remember !)
template <int R_A, int C_A, int C_B>
class mult_vv {
  int i, j, N, M, K;

 public:
  mult_vv(int i_, int j_, int N_, int M_, int K_)
      : i(i_), j(j_), N(N_), M(M_), K(K_) {}
  template <typename T>
  T operator()(Eigen::Matrix<T, -1, 1> x) const {
    using stan::math::multiply;
    Eigen::Matrix<T, R_A, C_A> A_c(N, M);
    Eigen::Matrix<T, C_A, C_B> B_c(M, K);
    int pos = 0;
    // traverse col-major
    for (int m = 0; m < M; ++m)
      for (int n = 0; n < N; ++n)
        A_c(n, m) = x(pos++);

    for (int k = 0; k < K; ++k)
      for (int m = 0; m < M; ++m)
        B_c(m, k) = x(pos++);

    Eigen::Matrix<T, R_A, C_B> AB_c = multiply(A_c, B_c);
    return AB_c(i, j);
  }
};

template <>
class mult_vv<1, -1, 1> {
  int N, M, K;

 public:
  mult_vv(int N_, int M_, int K_) : N(N_), M(M_), K(K_) {}
  template <typename T>
  T operator()(Eigen::Matrix<T, -1, 1> x) const {
    using stan::math::multiply;
    Eigen::Matrix<T, 1, -1> A_c(N, M);
    Eigen::Matrix<T, -1, 1> B_c(M, K);
    int pos = 0;
    // traverse col-major
    for (int m = 0; m < M; ++m)
      for (int n = 0; n < N; ++n)
        A_c(n, m) = x(pos++);

    for (int k = 0; k < K; ++k)
      for (int m = 0; m < M; ++m)
        B_c(m, k) = x(pos++);

    T AB_c = multiply(A_c, B_c);
    return AB_c;
  }
};

template <int R_A, int C_A, int C_B>
class mult_dv {
  int i, j, M, K;
  Eigen::Matrix<double, R_A, C_A> A_c;

 public:
  mult_dv(int i_, int j_, int M_, int K_, Eigen::Matrix<double, R_A, C_A> A_c_)
      : i(i_), j(j_), M(M_), K(K_), A_c(A_c_) {}
  template <typename T>
  T operator()(Eigen::Matrix<T, -1, 1> x) const {
    using stan::math::multiply;
    Eigen::Matrix<T, C_A, C_B> B_c(M, K);
    int pos = 0;
    // traverse col-major

    for (int k = 0; k < K; ++k)
      for (int m = 0; m < M; ++m)
        B_c(m, k) = x(pos++);

    Eigen::Matrix<T, R_A, C_B> AB_c = multiply(A_c, B_c);
    return AB_c(i, j);
  }
};

template <>
class mult_dv<1, -1, 1> {
  int M, K;
  Eigen::Matrix<double, 1, -1> A_c;

 public:
  mult_dv(int M_, int K_, Eigen::Matrix<double, 1, -1> A_c_)
      : M(M_), K(K_), A_c(A_c_) {}
  template <typename T>
  T operator()(Eigen::Matrix<T, -1, 1> x) const {
    using stan::math::multiply;
    Eigen::Matrix<T, -1, 1> B_c(M, K);
    int pos = 0;
    // traverse col-major

    for (int k = 0; k < K; ++k)
      for (int m = 0; m < M; ++m)
        B_c(m, k) = x(pos++);

    T AB_c = multiply(A_c, B_c);
    return AB_c;
  }
};

template <int R_A, int C_A, int C_B>
class mult_vd {
  int i, j, N, M;
  Eigen::Matrix<double, C_A, C_B> B_c;

 public:
  mult_vd(int i_, int j_, int N_, int M_, Eigen::Matrix<double, C_A, C_B> B_c_)
      : i(i_), j(j_), N(N_), M(M_), B_c(B_c_) {}
  template <typename T>
  T operator()(Eigen::Matrix<T, -1, 1> x) const {
    using stan::math::multiply;
    Eigen::Matrix<T, R_A, C_A> A_c(N, M);
    int pos = 0;
    // traverse col-major
    for (int m = 0; m < M; ++m)
      for (int n = 0; n < N; ++n)
        A_c(n, m) = x(pos++);

    Eigen::Matrix<T, -1, -1> AB_c = multiply(A_c, B_c);
    return AB_c(i, j);
  }
};

template <>
class mult_vd<1, -1, 1> {
  int N, M;
  Eigen::Matrix<double, -1, 1> B_c;

 public:
  mult_vd(int N_, int M_, Eigen::Matrix<double, -1, 1> B_c_)
      : N(N_), M(M_), B_c(B_c_) {}
  template <typename T>
  T operator()(Eigen::Matrix<T, -1, 1> x) const {
    using stan::math::multiply;
    Eigen::Matrix<T, 1, -1> A_c(N, M);
    int pos = 0;
    // traverse col-major
    for (int m = 0; m < M; ++m)
      for (int n = 0; n < N; ++n)
        A_c(n, m) = x(pos++);

    T AB_c = multiply(A_c, B_c);
    return AB_c;
  }
};

Eigen::Matrix<double, -1, 1> generate_inp(int N, int M, int K) {
  std::srand(123);
  int size_vec = N * M + M * K;
  Eigen::Matrix<double, -1, 1> vec
      = Eigen::Matrix<double, -1, 1>::Random(size_vec);
  return vec;
}

template <int R_A, int C_A, int C_B>
void pull_vals(int N, int M, int K, const Eigen::Matrix<double, -1, 1>& x,
               Eigen::Matrix<double, R_A, C_A>& A,
               Eigen::Matrix<double, C_A, C_B>& B) {
  A.resize(N, M);
  B.resize(M, K);
  int pos = 0;
  for (int m = 0; m < M; ++m)
    for (int n = 0; n < N; ++n)
      A(n, m) = x(pos++);

  for (int k = 0; k < K; ++k)
    for (int m = 0; m < M; ++m)
      B(m, k) = x(pos++);
}

TEST(AgradRevMatrix, multiply_scalar_scalar) {
  using stan::math::multiply;
  double d1, d2;
  AVAR v1, v2;

  d1 = 10;
  v1 = 10;
  d2 = -2;
  v2 = -2;

  EXPECT_FLOAT_EQ(-20.0, multiply(d1, d2));
  EXPECT_FLOAT_EQ(-20.0, multiply(d1, v2).val());
  EXPECT_FLOAT_EQ(-20.0, multiply(v1, d2).val());
  EXPECT_FLOAT_EQ(-20.0, multiply(v1, v2).val());

  EXPECT_FLOAT_EQ(6.0, multiply(AVAR(3), AVAR(2)).val());
  EXPECT_FLOAT_EQ(6.0, multiply(3.0, AVAR(2)).val());
  EXPECT_FLOAT_EQ(6.0, multiply(AVAR(3), 2.0).val());
}
TEST(AgradRevMatrix, multiply_vector_scalar) {
  using stan::math::vector_d;
  using stan::math::vector_v;

  vector_d d1(3);
  vector_v v1(3);
  double d2;
  AVAR v2;

  d1 << 100, 0, -3;
  v1 << 100, 0, -3;
  d2 = -2;
  v2 = -2;

  vector_v output;
  output = multiply(d1, v2);
  EXPECT_FLOAT_EQ(-200, output(0).val());
  EXPECT_FLOAT_EQ(0, output(1).val());
  EXPECT_FLOAT_EQ(6, output(2).val());

  output = multiply(v1, d2);
  EXPECT_FLOAT_EQ(-200, output(0).val());
  EXPECT_FLOAT_EQ(0, output(1).val());
  EXPECT_FLOAT_EQ(6, output(2).val());

  output = multiply(v1, v2);
  EXPECT_FLOAT_EQ(-200, output(0).val());
  EXPECT_FLOAT_EQ(0, output(1).val());
  EXPECT_FLOAT_EQ(6, output(2).val());
}
TEST(AgradRevMatrix, multiply_rowvector_scalar) {
  using stan::math::row_vector_d;
  using stan::math::row_vector_v;

  row_vector_d d1(3);
  row_vector_v v1(3);
  double d2;
  AVAR v2;

  d1 << 100, 0, -3;
  v1 << 100, 0, -3;
  d2 = -2;
  v2 = -2;

  row_vector_v output;
  output = multiply(d1, v2);
  EXPECT_FLOAT_EQ(-200, output(0).val());
  EXPECT_FLOAT_EQ(0, output(1).val());
  EXPECT_FLOAT_EQ(6, output(2).val());

  output = multiply(v1, d2);
  EXPECT_FLOAT_EQ(-200, output(0).val());
  EXPECT_FLOAT_EQ(0, output(1).val());
  EXPECT_FLOAT_EQ(6, output(2).val());

  output = multiply(v1, v2);
  EXPECT_FLOAT_EQ(-200, output(0).val());
  EXPECT_FLOAT_EQ(0, output(1).val());
  EXPECT_FLOAT_EQ(6, output(2).val());
}
TEST(AgradRevMatrix, multiply_matrix_scalar) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;

  matrix_d d1(2, 2);
  matrix_v v1(2, 2);
  double d2;
  AVAR v2;

  d1 << 100, 0, -3, 4;
  v1 << 100, 0, -3, 4;
  d2 = -2;
  v2 = -2;

  matrix_v output;
  output = multiply(d1, v2);
  EXPECT_FLOAT_EQ(-200, output(0, 0).val());
  EXPECT_FLOAT_EQ(0, output(0, 1).val());
  EXPECT_FLOAT_EQ(6, output(1, 0).val());
  EXPECT_FLOAT_EQ(-8, output(1, 1).val());

  output = multiply(v1, d2);
  EXPECT_FLOAT_EQ(-200, output(0, 0).val());
  EXPECT_FLOAT_EQ(0, output(0, 1).val());
  EXPECT_FLOAT_EQ(6, output(1, 0).val());
  EXPECT_FLOAT_EQ(-8, output(1, 1).val());

  output = multiply(v1, v2);
  EXPECT_FLOAT_EQ(-200, output(0, 0).val());
  EXPECT_FLOAT_EQ(0, output(0, 1).val());
  EXPECT_FLOAT_EQ(6, output(1, 0).val());
  EXPECT_FLOAT_EQ(-8, output(1, 1).val());
}
TEST(AgradRevMatrix, multiply_rowvector_vector) {
  using stan::math::row_vector_d;
  using stan::math::row_vector_v;
  using stan::math::vector_d;
  using stan::math::vector_v;

  row_vector_d d1(3);
  row_vector_v v1(3);
  vector_d d2(3);
  vector_v v2(3);

  d1 << 1, 3, -5;
  v1 << 1, 3, -5;
  d2 << 4, -2, -1;
  v2 << 4, -2, -1;

  EXPECT_FLOAT_EQ(3, multiply(v1, v2).val());
  EXPECT_FLOAT_EQ(3, multiply(v1, d2).val());
  EXPECT_FLOAT_EQ(3, multiply(d1, v2).val());

  d1.resize(1);
  v1.resize(1);
  EXPECT_THROW(multiply(v1, v2), std::invalid_argument);
  EXPECT_THROW(multiply(v1, d2), std::invalid_argument);
  EXPECT_THROW(multiply(d1, v2), std::invalid_argument);
}
TEST(AgradRevMatrix, multiply_vector_rowvector) {
  using stan::math::matrix_v;
  using stan::math::row_vector_d;
  using stan::math::row_vector_v;
  using stan::math::vector_d;
  using stan::math::vector_v;

  vector_d d1(3);
  vector_v v1(3);
  row_vector_d d2(3);
  row_vector_v v2(3);

  d1 << 1, 3, -5;
  v1 << 1, 3, -5;
  d2 << 4, -2, -1;
  v2 << 4, -2, -1;

  matrix_v output = multiply(v1, v2);
  EXPECT_EQ(3, output.rows());
  EXPECT_EQ(3, output.cols());
  EXPECT_FLOAT_EQ(4, output(0, 0).val());
  EXPECT_FLOAT_EQ(-2, output(0, 1).val());
  EXPECT_FLOAT_EQ(-1, output(0, 2).val());
  EXPECT_FLOAT_EQ(12, output(1, 0).val());
  EXPECT_FLOAT_EQ(-6, output(1, 1).val());
  EXPECT_FLOAT_EQ(-3, output(1, 2).val());
  EXPECT_FLOAT_EQ(-20, output(2, 0).val());
  EXPECT_FLOAT_EQ(10, output(2, 1).val());
  EXPECT_FLOAT_EQ(5, output(2, 2).val());

  output = multiply(v1, d2);
  EXPECT_EQ(3, output.rows());
  EXPECT_EQ(3, output.cols());
  EXPECT_FLOAT_EQ(4, output(0, 0).val());
  EXPECT_FLOAT_EQ(-2, output(0, 1).val());
  EXPECT_FLOAT_EQ(-1, output(0, 2).val());
  EXPECT_FLOAT_EQ(12, output(1, 0).val());
  EXPECT_FLOAT_EQ(-6, output(1, 1).val());
  EXPECT_FLOAT_EQ(-3, output(1, 2).val());
  EXPECT_FLOAT_EQ(-20, output(2, 0).val());
  EXPECT_FLOAT_EQ(10, output(2, 1).val());
  EXPECT_FLOAT_EQ(5, output(2, 2).val());

  output = multiply(d1, v2);
  EXPECT_EQ(3, output.rows());
  EXPECT_EQ(3, output.cols());
  EXPECT_FLOAT_EQ(4, output(0, 0).val());
  EXPECT_FLOAT_EQ(-2, output(0, 1).val());
  EXPECT_FLOAT_EQ(-1, output(0, 2).val());
  EXPECT_FLOAT_EQ(12, output(1, 0).val());
  EXPECT_FLOAT_EQ(-6, output(1, 1).val());
  EXPECT_FLOAT_EQ(-3, output(1, 2).val());
  EXPECT_FLOAT_EQ(-20, output(2, 0).val());
  EXPECT_FLOAT_EQ(10, output(2, 1).val());
  EXPECT_FLOAT_EQ(5, output(2, 2).val());
}
TEST(AgradRevMatrix, multiply_matrix_vector) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::vector_d;
  using stan::math::vector_v;

  matrix_d d1(3, 2);
  matrix_v v1(3, 2);
  vector_d d2(2);
  vector_v v2(2);

  d1 << 1, 3, -5, 4, -2, -1;
  v1 << 1, 3, -5, 4, -2, -1;
  d2 << -2, 4;
  v2 << -2, 4;

  vector_v output = multiply(v1, v2);
  EXPECT_EQ(3, output.size());
  EXPECT_FLOAT_EQ(10, output(0).val());
  EXPECT_FLOAT_EQ(26, output(1).val());
  EXPECT_FLOAT_EQ(0, output(2).val());

  output = multiply(v1, d2);
  EXPECT_EQ(3, output.size());
  EXPECT_FLOAT_EQ(10, output(0).val());
  EXPECT_FLOAT_EQ(26, output(1).val());
  EXPECT_FLOAT_EQ(0, output(2).val());

  output = multiply(d1, v2);
  EXPECT_EQ(3, output.size());
  EXPECT_FLOAT_EQ(10, output(0).val());
  EXPECT_FLOAT_EQ(26, output(1).val());
  EXPECT_FLOAT_EQ(0, output(2).val());
}
TEST(AgradRevMatrix, multiply_matrix_vector_exception) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::vector_d;
  using stan::math::vector_v;

  matrix_d d1(3, 2);
  matrix_v v1(3, 2);
  vector_d d2(4);
  vector_v v2(4);
  EXPECT_THROW(multiply(v1, v2), std::invalid_argument);
  EXPECT_THROW(multiply(v1, d2), std::invalid_argument);
  EXPECT_THROW(multiply(d1, v2), std::invalid_argument);
}
TEST(AgradRevMatrix, multiply_rowvector_matrix) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::row_vector_d;
  using stan::math::row_vector_v;
  using stan::math::vector_v;

  row_vector_d d1(3);
  row_vector_v v1(3);
  matrix_d d2(3, 2);
  matrix_v v2(3, 2);

  d1 << -2, 4, 1;
  v1 << -2, 4, 1;
  d2 << 1, 3, -5, 4, -2, -1;
  v2 << 1, 3, -5, 4, -2, -1;

  vector_v output = multiply(v1, v2);
  EXPECT_EQ(2, output.size());
  EXPECT_FLOAT_EQ(-24, output(0).val());
  EXPECT_FLOAT_EQ(9, output(1).val());

  output = multiply(v1, d2);
  EXPECT_EQ(2, output.size());
  EXPECT_FLOAT_EQ(-24, output(0).val());
  EXPECT_FLOAT_EQ(9, output(1).val());

  output = multiply(d1, v2);
  EXPECT_EQ(2, output.size());
  EXPECT_FLOAT_EQ(-24, output(0).val());
  EXPECT_FLOAT_EQ(9, output(1).val());
}
TEST(AgradRevMatrix, multiply_rowvector_matrix_exception) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::row_vector_d;
  using stan::math::row_vector_v;

  row_vector_d d1(4);
  row_vector_v v1(4);
  matrix_d d2(3, 2);
  matrix_v v2(3, 2);
  EXPECT_THROW(multiply(v1, v2), std::invalid_argument);
  EXPECT_THROW(multiply(v1, d2), std::invalid_argument);
  EXPECT_THROW(multiply(d1, v2), std::invalid_argument);
}
TEST(AgradRevMatrix, multiply_matrix_matrix) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;

  matrix_d d1(2, 3);
  matrix_v v1(2, 3);
  matrix_d d2(3, 2);
  matrix_v v2(3, 2);

  d1 << 9, 24, 3, 46, -9, -33;
  v1 << 9, 24, 3, 46, -9, -33;
  d2 << 1, 3, -5, 4, -2, -1;
  v2 << 1, 3, -5, 4, -2, -1;

  matrix_v output = multiply(v1, v2);
  EXPECT_EQ(2, output.rows());
  EXPECT_EQ(2, output.cols());
  EXPECT_FLOAT_EQ(-117, output(0, 0).val());
  EXPECT_FLOAT_EQ(120, output(0, 1).val());
  EXPECT_FLOAT_EQ(157, output(1, 0).val());
  EXPECT_FLOAT_EQ(135, output(1, 1).val());

  output = multiply(v1, d2);
  EXPECT_EQ(2, output.rows());
  EXPECT_EQ(2, output.cols());
  EXPECT_FLOAT_EQ(-117, output(0, 0).val());
  EXPECT_FLOAT_EQ(120, output(0, 1).val());
  EXPECT_FLOAT_EQ(157, output(1, 0).val());
  EXPECT_FLOAT_EQ(135, output(1, 1).val());

  output = multiply(d1, v2);
  EXPECT_EQ(2, output.rows());
  EXPECT_EQ(2, output.cols());
  EXPECT_FLOAT_EQ(-117, output(0, 0).val());
  EXPECT_FLOAT_EQ(120, output(0, 1).val());
  EXPECT_FLOAT_EQ(157, output(1, 0).val());
  EXPECT_FLOAT_EQ(135, output(1, 1).val());
}
TEST(AgradRevMatrix, multiply_matrix_matrix_exception) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;

  matrix_d d1(2, 2);
  matrix_v v1(2, 2);
  matrix_d d2(3, 2);
  matrix_v v2(3, 2);

  EXPECT_THROW(multiply(v1, v2), std::invalid_argument);
  EXPECT_THROW(multiply(v1, d2), std::invalid_argument);
  EXPECT_THROW(multiply(d1, v2), std::invalid_argument);
}
TEST(AgradRevMatrix, multiply_scalar_vector_cv) {
  using stan::math::multiply;
  using stan::math::vector_v;

  vector_v x(3);
  x << 1, 2, 3;
  AVEC x_ind = createAVEC(x(0), x(1), x(2));
  vector_v y = multiply(2.0, x);
  EXPECT_FLOAT_EQ(2.0, y(0).val());
  EXPECT_FLOAT_EQ(4.0, y(1).val());
  EXPECT_FLOAT_EQ(6.0, y(2).val());

  VEC g = cgradvec(y(0), x_ind);
  EXPECT_FLOAT_EQ(2.0, g[0]);
  EXPECT_FLOAT_EQ(0.0, g[1]);
  EXPECT_FLOAT_EQ(0.0, g[2]);
}
TEST(AgradRevMatrix, multiply_scalar_vector_vv) {
  using stan::math::multiply;
  using stan::math::vector_v;

  vector_v x(3);
  x << 1, 4, 9;
  AVAR two = 2.0;
  AVEC x_ind = createAVEC(x(0), x(1), x(2), two);
  vector_v y = multiply(two, x);
  EXPECT_FLOAT_EQ(2.0, y(0).val());
  EXPECT_FLOAT_EQ(8.0, y(1).val());
  EXPECT_FLOAT_EQ(18.0, y(2).val());

  VEC g = cgradvec(y(1), x_ind);
  EXPECT_FLOAT_EQ(0.0, g[0]);
  EXPECT_FLOAT_EQ(2.0, g[1]);
  EXPECT_FLOAT_EQ(0.0, g[2]);
  EXPECT_FLOAT_EQ(4.0, g[3]);
}
TEST(AgradRevMatrix, multiply_scalar_vector_vc) {
  using stan::math::multiply;
  using stan::math::vector_v;

  vector_v x(3);
  x << 1, 2, 3;
  AVAR two = 2.0;
  AVEC x_ind = createAVEC(two);
  vector_v y = multiply(two, x);
  EXPECT_FLOAT_EQ(2.0, y(0).val());
  EXPECT_FLOAT_EQ(4.0, y(1).val());
  EXPECT_FLOAT_EQ(6.0, y(2).val());

  VEC g = cgradvec(y(2), x_ind);
  EXPECT_FLOAT_EQ(3.0, g[0]);
}

TEST(AgradRevMatrix, multiply_scalar_row_vector_cv) {
  using stan::math::multiply;
  using stan::math::row_vector_v;

  row_vector_v x(3);
  x << 1, 2, 3;
  AVEC x_ind = createAVEC(x(0), x(1), x(2));
  row_vector_v y = multiply(2.0, x);
  EXPECT_FLOAT_EQ(2.0, y(0).val());
  EXPECT_FLOAT_EQ(4.0, y(1).val());
  EXPECT_FLOAT_EQ(6.0, y(2).val());

  VEC g = cgradvec(y(0), x_ind);
  EXPECT_FLOAT_EQ(2.0, g[0]);
  EXPECT_FLOAT_EQ(0.0, g[1]);
  EXPECT_FLOAT_EQ(0.0, g[2]);
}
TEST(AgradRevMatrix, multiply_scalar_row_vector_vv) {
  using stan::math::multiply;
  using stan::math::row_vector_v;

  row_vector_v x(3);
  x << 1, 4, 9;
  AVAR two = 2.0;
  AVEC x_ind = createAVEC(x(0), x(1), x(2), two);
  row_vector_v y = multiply(two, x);
  EXPECT_FLOAT_EQ(2.0, y(0).val());
  EXPECT_FLOAT_EQ(8.0, y(1).val());
  EXPECT_FLOAT_EQ(18.0, y(2).val());

  VEC g = cgradvec(y(1), x_ind);
  EXPECT_FLOAT_EQ(0.0, g[0]);
  EXPECT_FLOAT_EQ(2.0, g[1]);
  EXPECT_FLOAT_EQ(0.0, g[2]);
  EXPECT_FLOAT_EQ(4.0, g[3]);
}
TEST(AgradRevMatrix, multiply_scalar_row_vector_vc) {
  using stan::math::multiply;
  using stan::math::row_vector_v;

  row_vector_v x(3);
  x << 1, 2, 3;
  AVAR two = 2.0;
  AVEC x_ind = createAVEC(two);
  row_vector_v y = multiply(two, x);
  EXPECT_FLOAT_EQ(2.0, y(0).val());
  EXPECT_FLOAT_EQ(4.0, y(1).val());
  EXPECT_FLOAT_EQ(6.0, y(2).val());

  VEC g = cgradvec(y(2), x_ind);
  EXPECT_FLOAT_EQ(3.0, g[0]);
}

TEST(AgradRevMatrix, multiply_scalar_matrix_cv) {
  using stan::math::matrix_v;
  using stan::math::multiply;

  matrix_v x(2, 3);
  x << 1, 2, 3, 4, 5, 6;
  AVEC x_ind = createAVEC(x(0, 0), x(0, 1), x(0, 2), x(1, 0));
  matrix_v y = multiply(2.0, x);
  EXPECT_FLOAT_EQ(2.0, y(0, 0).val());
  EXPECT_FLOAT_EQ(4.0, y(0, 1).val());
  EXPECT_FLOAT_EQ(6.0, y(0, 2).val());

  VEC g = cgradvec(y(0, 0), x_ind);
  EXPECT_FLOAT_EQ(2.0, g[0]);
  EXPECT_FLOAT_EQ(0.0, g[1]);
  EXPECT_FLOAT_EQ(0.0, g[2]);
  EXPECT_FLOAT_EQ(0.0, g[3]);
}

TEST(AgradRevMatrix, multiply_scalar_matrix_vc) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::multiply;

  matrix_d x(2, 3);
  x << 1, 2, 3, 4, 5, 6;
  AVAR two = 2.0;
  AVEC x_ind = createAVEC(two);

  matrix_v y = multiply(two, x);
  EXPECT_FLOAT_EQ(2.0, y(0, 0).val());
  EXPECT_FLOAT_EQ(4.0, y(0, 1).val());
  EXPECT_FLOAT_EQ(6.0, y(0, 2).val());

  VEC g = cgradvec(y(1, 0), x_ind);
  EXPECT_FLOAT_EQ(4.0, g[0]);
}

TEST(AgradRevMatrix, multiply_vector_int) {
  // test namespace resolution
  using stan::math::multiply;
  using stan::math::multiply;
  using stan::math::vector_d;
  using stan::math::vector_v;

  vector_d dvec(3);
  dvec << 1, 2, 3;
  int a = 2;
  vector_d prod_vec = multiply(dvec, a);
  EXPECT_EQ(3, prod_vec.size());
  EXPECT_EQ(2.0, prod_vec[0]);
  EXPECT_EQ(4.0, prod_vec[1]);
  EXPECT_EQ(6.0, prod_vec[2]);
}

TEST(AgradRevMatrix, multiply_matrix_matrix_grad_fd) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  int N = 3;
  int M = 4;
  int K = 5;
  MatrixXd A;
  MatrixXd B;
  MatrixXd AB(N, K);
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_vv<-1, -1, -1> func(n, k, N, M, K);
      VectorXd grad_ad(N * M + M * K);
      VectorXd grad_fd(N * M + M * K);
      double val_ad;
      double val_fd;
      stan::math::gradient(func, test, val_ad, grad_ad);
      stan::math::finite_diff_gradient(func, test, val_fd, grad_fd);
      EXPECT_FLOAT_EQ(AB(n, k), val_ad);
      EXPECT_FLOAT_EQ(AB(n, k), val_fd);
      for (int i = 0; i < grad_ad.size(); ++i)
        EXPECT_NEAR(grad_ad(i), grad_fd(i), 1e-10);
    }
  }
}

TEST(AgradRevMatrix, multiply_matrix_matrix_grad_ex) {
  using Eigen::Infinity;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  int N = 3;
  int M = 4;
  int K = 5;
  MatrixXd A;
  MatrixXd B;
  MatrixXd AB(N, K);
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_vv<-1, -1, -1> func(n, k, N, M, K);
      VectorXd grad_ad(N * M + M * K);
      MatrixXd grad_A;
      MatrixXd grad_B;
      MatrixXd grad_A_ex;
      MatrixXd grad_B_ex;
      grad_A_ex.resize(N, M);
      grad_A_ex.setZero();
      grad_B_ex.resize(M, K);
      grad_B_ex.setZero();
      grad_A_ex.row(n) = B.col(k);
      grad_B_ex.col(k) = A.row(n);
      double val_ad;
      stan::math::gradient(func, test, val_ad, grad_ad);
      EXPECT_FLOAT_EQ(AB(n, k), val_ad);
      pull_vals(N, M, K, grad_ad, grad_A, grad_B);
      EXPECT_FLOAT_EQ((grad_A - grad_A_ex).lpNorm<Infinity>(), 0);
      EXPECT_FLOAT_EQ((grad_B - grad_B_ex).lpNorm<Infinity>(), 0);
    }
  }
}

TEST(AgradRevMatrix, multiply_matrix_vector_grad_fd) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  int N = 3;
  int M = 4;
  int K = 1;
  MatrixXd A;
  VectorXd B;
  VectorXd AB(N);
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_vv<-1, -1, 1> func(n, k, N, M, K);
      VectorXd grad_ad(N * M + M * K);
      VectorXd grad_fd(N * M + M * K);
      double val_ad;
      double val_fd;
      stan::math::gradient(func, test, val_ad, grad_ad);
      stan::math::finite_diff_gradient(func, test, val_fd, grad_fd);
      EXPECT_FLOAT_EQ(AB(n), val_ad);
      EXPECT_FLOAT_EQ(AB(n), val_fd);
      for (int i = 0; i < grad_ad.size(); ++i)
        EXPECT_NEAR(grad_ad(i), grad_fd(i), 1e-10);
    }
  }
}

TEST(AgradRevMatrix, multiply_matrix_vector_grad_ex) {
  using Eigen::Infinity;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  int N = 3;
  int M = 4;
  int K = 1;
  MatrixXd A;
  VectorXd B;
  VectorXd AB(N);
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_vv<-1, -1, 1> func(n, k, N, M, K);
      VectorXd grad_ad(N * M + M * K);
      MatrixXd grad_A;
      MatrixXd grad_B;
      MatrixXd grad_A_ex;
      MatrixXd grad_B_ex;
      grad_A_ex.resize(N, M);
      grad_A_ex.setZero();
      grad_B_ex.resize(M, K);
      grad_B_ex.setZero();
      grad_A_ex.row(n) = B.col(k);
      grad_B_ex.col(k) = A.row(n);
      double val_ad;
      stan::math::gradient(func, test, val_ad, grad_ad);
      EXPECT_FLOAT_EQ(AB(n), val_ad);
      pull_vals(N, M, K, grad_ad, grad_A, grad_B);
      EXPECT_FLOAT_EQ((grad_A - grad_A_ex).lpNorm<Infinity>(), 0);
      EXPECT_FLOAT_EQ((grad_B - grad_B_ex).lpNorm<Infinity>(), 0);
    }
  }
}

TEST(AgradRevMatrix, multiply_row_vector_matrix_grad_fd) {
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;

  int N = 1;
  int M = 4;
  int K = 5;
  RowVectorXd A;
  MatrixXd B;
  RowVectorXd AB(K);
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_vv<1, -1, -1> func(n, k, N, M, K);
      VectorXd grad_ad(N * M + M * K);
      VectorXd grad_fd(N * M + M * K);
      double val_ad;
      double val_fd;
      stan::math::gradient(func, test, val_ad, grad_ad);
      stan::math::finite_diff_gradient(func, test, val_fd, grad_fd);
      EXPECT_FLOAT_EQ(AB(k), val_ad);
      EXPECT_FLOAT_EQ(AB(k), val_fd);
      for (int i = 0; i < grad_ad.size(); ++i)
        EXPECT_NEAR(grad_ad(i), grad_fd(i), 1e-10);
    }
  }
}

TEST(AgradRevMatrix, multiply_row_vector_matrix_grad_ex) {
  using Eigen::Infinity;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;

  int N = 1;
  int M = 4;
  int K = 5;
  RowVectorXd A;
  MatrixXd B;
  RowVectorXd AB(K);
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_vv<1, -1, -1> func(n, k, N, M, K);
      VectorXd grad_ad(N * M + M * K);
      RowVectorXd grad_A;
      MatrixXd grad_B;
      RowVectorXd grad_A_ex;
      MatrixXd grad_B_ex;
      grad_A_ex.resize(M);
      grad_A_ex.setZero();
      grad_B_ex.resize(M, K);
      grad_B_ex.setZero();
      grad_A_ex.row(n) = B.col(k);
      grad_B_ex.col(k) = A.row(n);
      double val_ad;
      stan::math::gradient(func, test, val_ad, grad_ad);
      EXPECT_FLOAT_EQ(AB(k), val_ad);
      pull_vals(N, M, K, grad_ad, grad_A, grad_B);
      EXPECT_FLOAT_EQ((grad_A - grad_A_ex).lpNorm<Infinity>(), 0);
      EXPECT_FLOAT_EQ((grad_B - grad_B_ex).lpNorm<Infinity>(), 0);
    }
  }
}

TEST(AgradRevMatrix, multiply_row_vector_vector_grad_fd) {
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;

  int N = 1;
  int M = 4;
  int K = 1;
  RowVectorXd A;
  VectorXd B;
  double AB;
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_vv<1, -1, 1> func(N, M, K);
      VectorXd grad_ad(N * M + M * K);
      VectorXd grad_fd(N * M + M * K);
      double val_ad;
      double val_fd;
      stan::math::gradient(func, test, val_ad, grad_ad);
      stan::math::finite_diff_gradient(func, test, val_fd, grad_fd);
      EXPECT_FLOAT_EQ(AB, val_ad);
      EXPECT_FLOAT_EQ(AB, val_fd);
      for (int i = 0; i < grad_ad.size(); ++i)
        EXPECT_NEAR(grad_ad(i), grad_fd(i), 1e-10);
    }
  }
}

TEST(AgradRevMatrix, multiply_row_vector_vector_grad_ex) {
  using Eigen::Infinity;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;

  int N = 1;
  int M = 4;
  int K = 1;
  RowVectorXd A;
  VectorXd B;
  double AB;
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_vv<1, -1, 1> func(N, M, K);
      VectorXd grad_ad(N * M + M * K);
      RowVectorXd grad_A;
      VectorXd grad_B;
      RowVectorXd grad_A_ex;
      VectorXd grad_B_ex;
      grad_A_ex.resize(M);
      grad_A_ex.setZero();
      grad_B_ex.resize(M);
      grad_B_ex.setZero();
      grad_A_ex = B.transpose();
      grad_B_ex = A.transpose();
      double val_ad;
      stan::math::gradient(func, test, val_ad, grad_ad);
      EXPECT_FLOAT_EQ(AB, val_ad);
      pull_vals(N, M, K, grad_ad, grad_A, grad_B);
      EXPECT_FLOAT_EQ((grad_A - grad_A_ex).lpNorm<Infinity>(), 0);
      EXPECT_FLOAT_EQ((grad_B - grad_B_ex).lpNorm<Infinity>(), 0);
    }
  }
}

TEST(AgradRevMatrix, multiply_vector_row_vector_grad_fd) {
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;

  int N = 3;
  int M = 1;
  int K = 5;
  VectorXd A;
  RowVectorXd B;
  MatrixXd AB;
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_vv<-1, 1, -1> func(n, k, N, M, K);
      VectorXd grad_ad(N * M + M * K);
      VectorXd grad_fd(N * M + M * K);
      double val_ad;
      double val_fd;
      stan::math::gradient(func, test, val_ad, grad_ad);
      stan::math::finite_diff_gradient(func, test, val_fd, grad_fd);
      EXPECT_FLOAT_EQ(AB(n, k), val_ad);
      EXPECT_FLOAT_EQ(AB(n, k), val_fd);
      for (int i = 0; i < grad_ad.size(); ++i)
        EXPECT_NEAR(grad_ad(i), grad_fd(i), 1e-10);
    }
  }
}

TEST(AgradRevMatrix, multiply_vector_row_vector_grad_ex) {
  using Eigen::Infinity;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;

  int N = 3;
  int M = 1;
  int K = 5;
  VectorXd A;
  RowVectorXd B;
  MatrixXd AB;
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_vv<-1, 1, -1> func(n, k, N, M, K);
      VectorXd grad_ad(N * M + M * K);
      VectorXd grad_A;
      VectorXd grad_A_ex;
      RowVectorXd grad_B;
      RowVectorXd grad_B_ex;
      grad_A_ex.resize(N);
      grad_A_ex.setZero();
      grad_B_ex.resize(K);
      grad_B_ex.setZero();
      grad_A_ex.row(n) = B.col(k);
      grad_B_ex.col(k) = A.row(n);
      double val_ad;
      stan::math::gradient(func, test, val_ad, grad_ad);
      EXPECT_FLOAT_EQ(AB(n, k), val_ad);
      pull_vals(N, M, K, grad_ad, grad_A, grad_B);
      EXPECT_FLOAT_EQ((grad_A - grad_A_ex).lpNorm<Infinity>(), 0);
      EXPECT_FLOAT_EQ((grad_B - grad_B_ex).lpNorm<Infinity>(), 0);
    }
  }
}

TEST(AgradRevMatrix, multiply_matrix_matrix_grad_fd_dv) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  int N = 3;
  int M = 4;
  int K = 5;
  MatrixXd A;
  MatrixXd B;
  MatrixXd AB(N, K);
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_dv<-1, -1, -1> func(n, k, M, K, A);
      VectorXd grad_ad(M * K);
      VectorXd grad_fd(M * K);
      double val_ad;
      double val_fd;
      stan::math::gradient(func, test.tail(M * K), val_ad, grad_ad);
      stan::math::finite_diff_gradient(func, test.tail(M * K), val_fd, grad_fd);
      EXPECT_FLOAT_EQ(AB(n, k), val_ad);
      EXPECT_FLOAT_EQ(AB(n, k), val_fd);
      for (int i = 0; i < grad_ad.size(); ++i)
        EXPECT_NEAR(grad_ad(i), grad_fd(i), 1e-10);
    }
  }
}

TEST(AgradRevMatrix, multiply_matrix_matrix_grad_ex_dv) {
  using Eigen::Infinity;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  int N = 3;
  int M = 4;
  int K = 5;
  MatrixXd A;
  MatrixXd B;
  MatrixXd AB(N, K);
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_dv<-1, -1, -1> func(n, k, M, K, A);
      VectorXd grad_ad(M * K);
      MatrixXd grad_B;
      MatrixXd grad_A;
      MatrixXd grad_B_ex;
      grad_B_ex.resize(M, K);
      grad_B_ex.setZero();
      grad_B_ex.col(k) = A.row(n);
      double val_ad;
      stan::math::gradient(func, test.tail(M * K), val_ad, grad_ad);
      EXPECT_FLOAT_EQ(AB(n, k), val_ad);
      pull_vals(0, M, K, grad_ad, grad_A, grad_B);
      EXPECT_FLOAT_EQ((grad_B - grad_B_ex).lpNorm<Infinity>(), 0);
    }
  }
}

TEST(AgradRevMatrix, multiply_matrix_vector_grad_fd_dv) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  int N = 3;
  int M = 4;
  int K = 1;
  MatrixXd A;
  VectorXd B;
  VectorXd AB(N);
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_dv<-1, -1, 1> func(n, k, M, K, A);
      VectorXd grad_ad(M * K);
      VectorXd grad_fd(M * K);
      double val_ad;
      double val_fd;
      stan::math::gradient(func, test.tail(M * K), val_ad, grad_ad);
      stan::math::finite_diff_gradient(func, test.tail(M * K), val_fd, grad_fd);
      EXPECT_FLOAT_EQ(AB(n), val_ad);
      EXPECT_FLOAT_EQ(AB(n), val_fd);
      for (int i = 0; i < grad_ad.size(); ++i)
        EXPECT_NEAR(grad_ad(i), grad_fd(i), 1e-10);
    }
  }
}

TEST(AgradRevMatrix, multiply_matrix_vector_grad_ex_dv) {
  using Eigen::Infinity;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  int N = 3;
  int M = 4;
  int K = 1;
  MatrixXd A;
  VectorXd B;
  VectorXd AB(N);
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_dv<-1, -1, 1> func(n, k, M, K, A);
      VectorXd grad_ad(M * K);
      VectorXd grad_B;
      MatrixXd grad_A;
      VectorXd grad_B_ex;
      grad_B_ex.resize(M);
      grad_B_ex.setZero();
      grad_B_ex = A.row(n);
      double val_ad;
      stan::math::gradient(func, test.tail(M * K), val_ad, grad_ad);
      EXPECT_FLOAT_EQ(AB(n), val_ad);
      pull_vals(0, M, K, grad_ad, grad_A, grad_B);
      EXPECT_FLOAT_EQ((grad_B - grad_B_ex).lpNorm<Infinity>(), 0);
    }
  }
}

TEST(AgradRevMatrix, multiply_row_vector_matrix_grad_fd_dv) {
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;

  int N = 1;
  int M = 4;
  int K = 5;
  RowVectorXd A;
  MatrixXd B;
  RowVectorXd AB(K);
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_dv<1, -1, -1> func(n, k, M, K, A);
      VectorXd grad_ad(M * K);
      VectorXd grad_fd(M * K);
      double val_ad;
      double val_fd;
      stan::math::gradient(func, test.tail(M * K), val_ad, grad_ad);
      stan::math::finite_diff_gradient(func, test.tail(M * K), val_fd, grad_fd);
      EXPECT_FLOAT_EQ(AB(k), val_ad);
      EXPECT_FLOAT_EQ(AB(k), val_fd);
      for (int i = 0; i < grad_ad.size(); ++i)
        EXPECT_NEAR(grad_ad(i), grad_fd(i), 1e-10);
    }
  }
}

TEST(AgradRevMatrix, multiply_row_vector_matrix_grad_ex_dv) {
  using Eigen::Infinity;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;

  int N = 1;
  int M = 4;
  int K = 5;
  RowVectorXd A;
  MatrixXd B;
  RowVectorXd AB(N, K);
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_dv<1, -1, -1> func(n, k, M, K, A);
      VectorXd grad_ad(M * K);
      MatrixXd grad_B;
      MatrixXd grad_A;
      MatrixXd grad_B_ex;
      grad_B_ex.resize(M, K);
      grad_B_ex.setZero();
      grad_B_ex.col(k) = A.row(n);
      double val_ad;
      stan::math::gradient(func, test.tail(M * K), val_ad, grad_ad);
      EXPECT_FLOAT_EQ(AB(k), val_ad);
      pull_vals(0, M, K, grad_ad, grad_A, grad_B);
      EXPECT_FLOAT_EQ((grad_B - grad_B_ex).lpNorm<Infinity>(), 0);
    }
  }
}

TEST(AgradRevMatrix, multiply_row_vector_vector_grad_fd_dv) {
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;

  int N = 1;
  int M = 4;
  int K = 1;
  RowVectorXd A;
  VectorXd B;
  double AB;
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_dv<1, -1, 1> func(M, K, A);
      VectorXd grad_ad(M * K);
      VectorXd grad_fd(M * K);
      double val_ad;
      double val_fd;
      stan::math::gradient(func, test.tail(M * K), val_ad, grad_ad);
      stan::math::finite_diff_gradient(func, test.tail(M * K), val_fd, grad_fd);
      EXPECT_FLOAT_EQ(AB, val_ad);
      EXPECT_FLOAT_EQ(AB, val_fd);
      for (int i = 0; i < grad_ad.size(); ++i)
        EXPECT_NEAR(grad_ad(i), grad_fd(i), 1e-10);
    }
  }
}

TEST(AgradRevMatrix, multiply_row_vector_vector_grad_ex_dv) {
  using Eigen::Infinity;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;

  int N = 1;
  int M = 4;
  int K = 1;
  RowVectorXd A;
  VectorXd B;
  double AB;
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_dv<1, -1, 1> func(M, K, A);
      VectorXd grad_ad(M * K);
      VectorXd grad_B;
      MatrixXd grad_A;
      VectorXd grad_B_ex;
      grad_B_ex.resize(M, K);
      grad_B_ex.setZero();
      grad_B_ex = A;
      double val_ad;
      stan::math::gradient(func, test.tail(M * K), val_ad, grad_ad);
      EXPECT_FLOAT_EQ(AB, val_ad);
      pull_vals(0, M, K, grad_ad, grad_A, grad_B);
      EXPECT_FLOAT_EQ((grad_B - grad_B_ex).lpNorm<Infinity>(), 0);
    }
  }
}

TEST(AgradRevMatrix, multiply_vector_row_vector_grad_fd_dv) {
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;

  int N = 3;
  int M = 1;
  int K = 5;
  VectorXd A;
  RowVectorXd B;
  MatrixXd AB;
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_dv<-1, 1, -1> func(n, k, M, K, A);
      VectorXd grad_ad(M * K);
      VectorXd grad_fd(M * K);
      double val_ad;
      double val_fd;
      stan::math::gradient(func, test.tail(M * K), val_ad, grad_ad);
      stan::math::finite_diff_gradient(func, test.tail(M * K), val_fd, grad_fd);
      EXPECT_FLOAT_EQ(AB(n, k), val_ad);
      EXPECT_FLOAT_EQ(AB(n, k), val_fd);
      for (int i = 0; i < grad_ad.size(); ++i)
        EXPECT_NEAR(grad_ad(i), grad_fd(i), 1e-10);
    }
  }
}

TEST(AgradRevMatrix, multiply_vector_row_vector_grad_ex_dv) {
  using Eigen::Infinity;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;

  int N = 3;
  int M = 1;
  int K = 5;
  VectorXd A;
  RowVectorXd B;
  MatrixXd AB;
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_dv<-1, 1, -1> func(n, k, M, K, A);
      VectorXd grad_ad(M * K);
      VectorXd grad_A;
      RowVectorXd grad_B;
      RowVectorXd grad_B_ex;
      grad_B_ex.resize(K);
      grad_B_ex.setZero();
      grad_B_ex.col(k) = A.row(n);
      double val_ad;
      stan::math::gradient(func, test.tail(M * K), val_ad, grad_ad);
      EXPECT_FLOAT_EQ(AB(n, k), val_ad);
      pull_vals(0, M, K, grad_ad, grad_A, grad_B);
      EXPECT_FLOAT_EQ((grad_B - grad_B_ex).lpNorm<Infinity>(), 0);
    }
  }
}

TEST(AgradRevMatrix, multiply_matrix_matrix_grad_fd_vd) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  int N = 3;
  int M = 4;
  int K = 5;
  MatrixXd A;
  MatrixXd B;
  MatrixXd AB(N, K);
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_vd<-1, -1, -1> func(n, k, N, M, B);
      VectorXd grad_ad(N * M);
      VectorXd grad_fd(N * M);
      double val_ad;
      double val_fd;
      stan::math::gradient(func, test.head(N * M), val_ad, grad_ad);
      stan::math::finite_diff_gradient(func, test.head(N * M), val_fd, grad_fd);
      EXPECT_FLOAT_EQ(AB(n, k), val_ad);
      EXPECT_FLOAT_EQ(AB(n, k), val_fd);
      for (int i = 0; i < grad_ad.size(); ++i)
        EXPECT_NEAR(grad_ad(i), grad_fd(i), 1e-10);
    }
  }
}

TEST(AgradRevMatrix, multiply_matrix_matrix_grad_ex_vd) {
  using Eigen::Infinity;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  int N = 3;
  int M = 4;
  int K = 5;
  MatrixXd A;
  MatrixXd B;
  MatrixXd AB(N, K);
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_vd<-1, -1, -1> func(n, k, N, M, B);
      VectorXd grad_ad(N * M);
      MatrixXd grad_B;
      MatrixXd grad_A;
      MatrixXd grad_A_ex;
      grad_A_ex.resize(N, M);
      grad_A_ex.setZero();
      grad_A_ex.row(n) = B.col(k);
      double val_ad;
      stan::math::gradient(func, test.head(N * M), val_ad, grad_ad);
      EXPECT_FLOAT_EQ(AB(n, k), val_ad);
      pull_vals(N, M, 0, grad_ad, grad_A, grad_B);
      EXPECT_FLOAT_EQ((grad_A - grad_A_ex).lpNorm<Infinity>(), 0);
    }
  }
}

TEST(AgradRevMatrix, multiply_matrix_vector_grad_fd_vd) {
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  int N = 3;
  int M = 4;
  int K = 1;
  MatrixXd A;
  VectorXd B;
  VectorXd AB(N);
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_vd<-1, -1, 1> func(n, k, N, M, B);
      VectorXd grad_ad(N * M);
      VectorXd grad_fd(N * M);
      double val_ad;
      double val_fd;
      stan::math::gradient(func, test.head(M * N), val_ad, grad_ad);
      stan::math::finite_diff_gradient(func, test.head(N * M), val_fd, grad_fd);
      EXPECT_FLOAT_EQ(AB(n), val_ad);
      EXPECT_FLOAT_EQ(AB(n), val_fd);
      for (int i = 0; i < grad_ad.size(); ++i)
        EXPECT_NEAR(grad_ad(i), grad_fd(i), 1e-10);
    }
  }
}

TEST(AgradRevMatrix, multiply_matrix_vector_grad_ex_vd) {
  using Eigen::Infinity;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  int N = 3;
  int M = 4;
  int K = 1;
  MatrixXd A;
  MatrixXd B;
  VectorXd AB(N);
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_vd<-1, -1, 1> func(n, k, N, M, B);
      VectorXd grad_ad(N * M);
      MatrixXd grad_B;
      MatrixXd grad_A;
      MatrixXd grad_A_ex;
      grad_A_ex.resize(N, M);
      grad_A_ex.setZero();
      grad_A_ex.row(n) = B.col(k);
      double val_ad;
      stan::math::gradient(func, test.head(N * M), val_ad, grad_ad);
      EXPECT_FLOAT_EQ(AB(n), val_ad);
      pull_vals(N, M, 0, grad_ad, grad_A, grad_B);
      EXPECT_FLOAT_EQ((grad_A - grad_A_ex).lpNorm<Infinity>(), 0);
    }
  }
}

TEST(AgradRevMatrix, multiply_row_vector_matrix_grad_fd_vd) {
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;

  int N = 1;
  int M = 4;
  int K = 5;
  RowVectorXd A;
  MatrixXd B;
  RowVectorXd AB(K);
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_vd<1, -1, -1> func(n, k, N, M, B);
      VectorXd grad_ad(N * M);
      VectorXd grad_fd(N * M);
      double val_ad;
      double val_fd;
      stan::math::gradient(func, test.head(N * M), val_ad, grad_ad);
      stan::math::finite_diff_gradient(func, test.head(N * M), val_fd, grad_fd);
      EXPECT_FLOAT_EQ(AB(k), val_ad);
      EXPECT_FLOAT_EQ(AB(k), val_fd);
      for (int i = 0; i < grad_ad.size(); ++i)
        EXPECT_NEAR(grad_ad(i), grad_fd(i), 1e-10);
    }
  }
}

TEST(AgradRevMatrix, multiply_row_vector_matrix_grad_ex_vd) {
  using Eigen::Infinity;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;

  int N = 1;
  int M = 4;
  int K = 5;
  RowVectorXd A;
  MatrixXd B;
  RowVectorXd AB(K);
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_vd<1, -1, -1> func(n, k, N, M, B);
      VectorXd grad_ad(N * M);
      MatrixXd grad_B;
      MatrixXd grad_A;
      RowVectorXd grad_A_ex;
      grad_A_ex.resize(M);
      grad_A_ex.setZero();
      grad_A_ex = B.col(k);
      double val_ad;
      stan::math::gradient(func, test.head(N * M), val_ad, grad_ad);
      EXPECT_FLOAT_EQ(AB(k), val_ad);
      pull_vals(N, M, 0, grad_ad, grad_A, grad_B);
      EXPECT_FLOAT_EQ((grad_A - grad_A_ex).lpNorm<Infinity>(), 0);
    }
  }
}

TEST(AgradRevMatrix, multiply_row_vector_vector_grad_fd_vd) {
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;

  int N = 1;
  int M = 4;
  int K = 1;
  RowVectorXd A;
  VectorXd B;
  double AB;
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_vd<1, -1, 1> func(N, M, B);
      VectorXd grad_ad(N * M);
      VectorXd grad_fd(N * M);
      double val_ad;
      double val_fd;
      stan::math::gradient(func, test.head(N * M), val_ad, grad_ad);
      stan::math::finite_diff_gradient(func, test.head(N * M), val_fd, grad_fd);
      EXPECT_FLOAT_EQ(AB, val_ad);
      EXPECT_FLOAT_EQ(AB, val_fd);
      for (int i = 0; i < grad_ad.size(); ++i)
        EXPECT_NEAR(grad_ad(i), grad_fd(i), 1e-10);
    }
  }
}

TEST(AgradRevMatrix, multiply_row_vector_vector_grad_ex_vd) {
  using Eigen::Infinity;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;

  int N = 1;
  int M = 4;
  int K = 1;
  RowVectorXd A;
  VectorXd B;
  double AB;
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_vd<1, -1, 1> func(N, M, B);
      VectorXd grad_ad(N * M);
      MatrixXd grad_B;
      RowVectorXd grad_A;
      RowVectorXd grad_A_ex;
      grad_A_ex.resize(M);
      grad_A_ex.setZero();
      grad_A_ex = B;
      double val_ad;
      stan::math::gradient(func, test.head(N * M), val_ad, grad_ad);
      EXPECT_FLOAT_EQ(AB, val_ad);
      pull_vals(N, M, 0, grad_ad, grad_A, grad_B);
      EXPECT_FLOAT_EQ((grad_A - grad_A_ex).lpNorm<Infinity>(), 0);
    }
  }
}

TEST(AgradRevMatrix, multiply_vector_row_vector_grad_fd_vd) {
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;

  int N = 3;
  int M = 1;
  int K = 5;
  VectorXd A;
  RowVectorXd B;
  MatrixXd AB;
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_vd<-1, 1, -1> func(n, k, N, M, B);
      VectorXd grad_ad(N * M);
      VectorXd grad_fd(N * M);
      double val_ad;
      double val_fd;
      stan::math::gradient(func, test.head(M * N), val_ad, grad_ad);
      stan::math::finite_diff_gradient(func, test.head(N * M), val_fd, grad_fd);
      EXPECT_FLOAT_EQ(AB(n, k), val_ad);
      EXPECT_FLOAT_EQ(AB(n, k), val_fd);
      for (int i = 0; i < grad_ad.size(); ++i)
        EXPECT_NEAR(grad_ad(i), grad_fd(i), 1e-10);
    }
  }
}

TEST(AgradRevMatrix, multiply_vector_row_vector_grad_ex_vd) {
  using Eigen::Infinity;
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;

  int N = 3;
  int M = 1;
  int K = 5;
  VectorXd A;
  RowVectorXd B;
  MatrixXd AB;
  VectorXd test = generate_inp(N, M, K);
  pull_vals(N, M, K, test, A, B);
  AB = A * B;
  for (int n = 0; n < N; ++n) {
    for (int k = 0; k < K; ++k) {
      mult_vd<-1, 1, -1> func(n, k, N, M, B);
      VectorXd grad_ad(N * M);
      RowVectorXd grad_B;
      VectorXd grad_A;
      VectorXd grad_A_ex;
      grad_A_ex.resize(N);
      grad_A_ex.setZero();
      grad_A_ex.row(n) = B.col(k);
      double val_ad;
      stan::math::gradient(func, test.head(N * M), val_ad, grad_ad);
      EXPECT_FLOAT_EQ(AB(n, k), val_ad);
      pull_vals(N, M, 0, grad_ad, grad_A, grad_B);
      EXPECT_FLOAT_EQ((grad_A - grad_A_ex).lpNorm<Infinity>(), 0);
    }
  }
}
TEST(AgradRevMatrix, check_varis_on_stack) {
  using stan::math::value_of;
  stan::math::matrix_v m(3, 3);
  m << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  stan::math::vector_v v(3);
  v << 10, 20, 30;
  stan::math::row_vector_v rv(3);
  rv << 100, 200, 300;
  stan::math::var s = 1;

  test::check_varis_on_stack(stan::math::multiply(m, m));
  test::check_varis_on_stack(stan::math::multiply(m, value_of(m)));
  test::check_varis_on_stack(stan::math::multiply(value_of(m), m));

  test::check_varis_on_stack(stan::math::multiply(m, v));
  test::check_varis_on_stack(stan::math::multiply(m, value_of(v)));
  test::check_varis_on_stack(stan::math::multiply(value_of(m), v));

  test::check_varis_on_stack(stan::math::multiply(rv, m));
  test::check_varis_on_stack(stan::math::multiply(rv, value_of(m)));
  test::check_varis_on_stack(stan::math::multiply(value_of(rv), m));

  test::check_varis_on_stack(stan::math::multiply(rv, v));
  test::check_varis_on_stack(stan::math::multiply(rv, value_of(v)));
  test::check_varis_on_stack(stan::math::multiply(value_of(rv), v));

  test::check_varis_on_stack(stan::math::multiply(s, m));
  test::check_varis_on_stack(stan::math::multiply(s, value_of(m)));
  test::check_varis_on_stack(stan::math::multiply(value_of(s), m));

  test::check_varis_on_stack(stan::math::multiply(s, rv));
  test::check_varis_on_stack(stan::math::multiply(s, value_of(rv)));
  test::check_varis_on_stack(stan::math::multiply(value_of(s), rv));

  test::check_varis_on_stack(stan::math::multiply(s, v));
  test::check_varis_on_stack(stan::math::multiply(s, value_of(v)));
  test::check_varis_on_stack(stan::math::multiply(value_of(s), v));

  test::check_varis_on_stack(stan::math::multiply(m, s));
  test::check_varis_on_stack(stan::math::multiply(m, value_of(s)));
  test::check_varis_on_stack(stan::math::multiply(value_of(m), s));

  test::check_varis_on_stack(stan::math::multiply(rv, s));
  test::check_varis_on_stack(stan::math::multiply(rv, value_of(s)));
  test::check_varis_on_stack(stan::math::multiply(value_of(rv), s));

  test::check_varis_on_stack(stan::math::multiply(v, s));
  test::check_varis_on_stack(stan::math::multiply(v, value_of(s)));
  test::check_varis_on_stack(stan::math::multiply(value_of(v), s));

  test::check_varis_on_stack(stan::math::multiply(s, s));
  test::check_varis_on_stack(stan::math::multiply(s, value_of(s)));
  test::check_varis_on_stack(stan::math::multiply(value_of(s), s));
}
