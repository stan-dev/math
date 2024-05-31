#include <stan/math.hpp>
#include <test/unit/util.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>

TEST_F(AgradRev, arena_matrix_matrix_test) {
  using Eigen::MatrixXd;
  using stan::math::arena_matrix;

  // construction
  arena_matrix<MatrixXd> a;
  arena_matrix<MatrixXd> a2;
  arena_matrix<MatrixXd> a3;
  arena_matrix<MatrixXd> b(3, 2);
  arena_matrix<MatrixXd> b2(4, 5);
  arena_matrix<MatrixXd> c(MatrixXd::Ones(3, 2));
  arena_matrix<MatrixXd> d(c);
  arena_matrix<MatrixXd> e(2 * d);

  // assignment
  a = c;
  a2 = std::move(d);
  a3 = 2 * a;
  b = d;
  b2 = std::move(c);
  e = e + a;
  a = MatrixXd::Ones(3, 2);

  EXPECT_MATRIX_EQ(a + a2 + a3 + b + b2 + e, MatrixXd::Ones(3, 2) * 9);
}

TEST_F(AgradRev, arena_matrix_vector_test) {
  using Eigen::VectorXd;
  using stan::math::arena_matrix;

  // construction
  arena_matrix<VectorXd> a;
  arena_matrix<VectorXd> a2;
  arena_matrix<VectorXd> a3;
  arena_matrix<VectorXd> b(3);
  arena_matrix<VectorXd> b2(4);
  arena_matrix<VectorXd> c(VectorXd::Ones(3));
  arena_matrix<VectorXd> d(c);
  arena_matrix<VectorXd> e(2 * d);

  // assignment
  a = c;
  a2 = std::move(d);
  a3 = 2 * a;
  b = d;
  b2 = std::move(c);
  e = e + a;
  a = VectorXd::Ones(3);

  EXPECT_MATRIX_EQ(a + a2 + a3 + b + b2 + e, VectorXd::Ones(3) * 9);
}

TEST_F(AgradRev, arena_matrix_row_vector_test) {
  using Eigen::RowVectorXd;
  using stan::math::arena_matrix;

  // construction
  arena_matrix<RowVectorXd> a;
  arena_matrix<RowVectorXd> a2;
  arena_matrix<RowVectorXd> a3;
  arena_matrix<RowVectorXd> b(3);
  arena_matrix<RowVectorXd> b2(4);
  arena_matrix<RowVectorXd> c(RowVectorXd::Ones(3));
  arena_matrix<RowVectorXd> d(c);
  arena_matrix<RowVectorXd> e(2 * d);

  // assignment
  a = c;
  a2 = std::move(d);
  a3 = 2 * a;
  b = d;
  b2 = std::move(c);
  e = e + a;
  a = RowVectorXd::Ones(3);

  EXPECT_MATRIX_EQ(a + a2 + a3 + b + b2 + e, RowVectorXd::Ones(3) * 9);
}

TEST_F(AgradRev, arena_matrix_transpose_test) {
  using stan::math::arena_matrix;

  Eigen::VectorXd c = Eigen::VectorXd::Random(3);
  Eigen::RowVectorXd d = Eigen::RowVectorXd::Random(3);

  // The VectorXd/RowVectorXd mixup here is on purpose to test a vector can
  // initialize a row vector and visa versa
  arena_matrix<Eigen::VectorXd> a = d;
  arena_matrix<Eigen::RowVectorXd> b = c;
  EXPECT_MATRIX_EQ(a, d.transpose());
  EXPECT_MATRIX_EQ(b, c.transpose());

  EXPECT_NO_THROW(a = b);
  EXPECT_MATRIX_EQ(a, b.transpose());
  a = Eigen::VectorXd::Random(3);
  EXPECT_NO_THROW(b = a);
  EXPECT_MATRIX_EQ(a, b.transpose());
}

template <typename T_map, typename T_mat>
inline void expect_sparse_matrix_equal(T_map&& map_x, T_mat&& mat_x) {
  using inner_iterator = typename std::decay_t<T_mat>::InnerIterator;
  using map_inner_iterator = typename std::decay_t<T_map>::InnerIterator;
  EXPECT_EQ(map_x.nonZeros(), mat_x.nonZeros());
  for (int k = 0; k < map_x.outerSize(); ++k) {
    inner_iterator it(mat_x, k);
    map_inner_iterator iz(map_x, k);
    for (; it && iz; ++it, ++iz) {
      EXPECT_FLOAT_EQ(iz.value(), it.value());
    }
  }
}

template <typename T_map>
inline void expect_sparse_dense_matrix_equal(
    T_map&& map_x, const Eigen::Matrix<double, -1, -1>& mat_x) {
  using map_inner_iterator = typename std::decay_t<T_map>::InnerIterator;
  Eigen::Index nnz_m = map_x.nonZeros();
  for (int k = 0; k < map_x.outerSize(); ++k) {
    for (map_inner_iterator iz(map_x, k); iz; ++iz) {
      EXPECT_FLOAT_EQ(iz.value(), mat_x(iz.row(), iz.col()));
    }
  }
}

TEST_F(AgradRev, arena_sparse_matrix_constructors) {
  using eig_mat = Eigen::SparseMatrix<double>;
  using stan::math::arena_matrix;
  using arena_mat = arena_matrix<eig_mat>;
  using inner_iterator = typename eig_mat::InnerIterator;
  using map_inner_iterator = typename Eigen::Map<eig_mat>::InnerIterator;
  using stan::test::make_sparse_matrix_random;
  eig_mat A = make_sparse_matrix_random(10, 10);
  // Testing each constructor
  arena_mat empty_m{};
  arena_mat nocopy(A.rows(), A.cols(), A.nonZeros(), A.outerIndexPtr(),
                   A.innerIndexPtr(), A.valuePtr(), A.innerNonZeroPtr());
  Eigen::Map<eig_mat> sparse_map(A.rows(), A.cols(), A.nonZeros(),
                                 A.outerIndexPtr(), A.innerIndexPtr(),
                                 A.valuePtr(), A.innerNonZeroPtr());
  arena_mat sparse_map_constructor(sparse_map);
  arena_mat A_m(A);
  expect_sparse_matrix_equal(A_m, A);
  arena_mat arena_mat_constructor(A_m);
  expect_sparse_matrix_equal(arena_mat_constructor, A);
  eig_mat mul_A = A * A;
  arena_mat mul_A_m(A * A_m);
  expect_sparse_matrix_equal(mul_A_m, mul_A);
  eig_mat mul_B = mul_A * A;
  arena_mat mul_B_m(mul_A_m * A_m);
  expect_sparse_matrix_equal(mul_B_m, mul_B);
}

TEST_F(AgradRev, arena_sparse_matrix_equals) {
  using eig_mat = Eigen::SparseMatrix<double>;
  using stan::math::arena_matrix;
  using arena_mat = arena_matrix<eig_mat>;
  using inner_iterator = typename eig_mat::InnerIterator;
  using map_inner_iterator = typename Eigen::Map<eig_mat>::InnerIterator;
  using stan::test::make_sparse_matrix_random;
  eig_mat A = make_sparse_matrix_random(10, 10);
  // Testing each constructor
  arena_mat empty_m = arena_mat{};
  empty_m = arena_mat{};
  arena_mat nocopy
      = arena_mat(A.rows(), A.cols(), A.nonZeros(), A.outerIndexPtr(),
                  A.innerIndexPtr(), A.valuePtr(), A.innerNonZeroPtr());
  nocopy = arena_mat(A.rows(), A.cols(), A.nonZeros(), A.outerIndexPtr(),
                     A.innerIndexPtr(), A.valuePtr(), A.innerNonZeroPtr());
  Eigen::Map<eig_mat> sparse_map(A.rows(), A.cols(), A.nonZeros(),
                                 A.outerIndexPtr(), A.innerIndexPtr(),
                                 A.valuePtr(), A.innerNonZeroPtr());
  arena_mat sparse_map_constructor = sparse_map;
  sparse_map_constructor = sparse_map;
  arena_mat A_m = A;
  A_m = A;
  expect_sparse_matrix_equal(A_m, A);
  arena_mat arena_mat_constructor = A_m;
  arena_mat_constructor = A_m;
  expect_sparse_matrix_equal(arena_mat_constructor, A);
  eig_mat mul_A = A * A;
  arena_mat mul_A_m = A * A_m;
  mul_A_m = A * A_m;
  expect_sparse_matrix_equal(mul_A_m, mul_A);
  eig_mat mul_B = mul_A * A;
  arena_mat mul_B_m = mul_A_m * A_m;
  mul_B_m = mul_A_m * A_m;
  expect_sparse_matrix_equal(mul_B_m, mul_B);
}

TEST_F(AgradRev, arena_sparse_matrix_inplace_ops) {
  using eig_mat = Eigen::SparseMatrix<double>;
  using stan::math::arena_matrix;
  using arena_mat = arena_matrix<eig_mat>;
  using inner_iterator = typename eig_mat::InnerIterator;
  using map_inner_iterator = typename Eigen::Map<eig_mat>::InnerIterator;
  using stan::test::make_sparse_matrix_random;
  eig_mat A = make_sparse_matrix_random(10, 10);
  arena_mat A_m = A;
  A_m += A;
  expect_sparse_matrix_equal(A_m, eig_mat(A + A));

  A_m = A;
  A_m -= A;
  expect_sparse_matrix_equal(A_m, eig_mat(A - A));

  arena_mat B_m = A;
  A_m = A;
  A_m += B_m;
  expect_sparse_matrix_equal(A_m, eig_mat(A + A));

  A_m = A;
  A_m -= B_m;
  expect_sparse_matrix_equal(A_m, eig_mat(A - A));

  Eigen::Matrix<double, -1, -1> B
      = Eigen::Matrix<double, -1, -1>::Random(10, 10);
  A_m = A;
  A_m += B;
  expect_sparse_dense_matrix_equal(A_m, Eigen::Matrix<double, -1, -1>(A + B));

  A_m = A;
  A_m -= B;
  expect_sparse_dense_matrix_equal(A_m, Eigen::Matrix<double, -1, -1>(A - B));

  double c = 10;
  A_m = A;
  A_m += c;
  eig_mat C = A;
  for (int k = 0; k < C.outerSize(); ++k) {
    inner_iterator it(C, k);
    for (; static_cast<bool>(it); ++it) {
      it.valueRef() += c;
    }
  }
  expect_sparse_dense_matrix_equal(A_m, C);
}
