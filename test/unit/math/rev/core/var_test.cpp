#include <stan/math.hpp>
#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/math/rev/core/gradable.hpp>
#include <gtest/gtest.h>
#include <string>
#include <vector>

struct AgradRev : public testing::Test {
  void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
  }
};

namespace stan {
namespace test {
template <typename T, typename S>
void ctor_overloads_float_impl() {
  using stan::math::var_value;
  using stan::math::vari_value;
  using stan::math::test::type_name;
  // standard constructor
  EXPECT_FLOAT_EQ(3.7, var_value<T>(3.7).val())
      << "Failed For T: " << type_name<T>() << "\n";
  // make sure copy ctor is used rather than casting vari* to unsigned int
  EXPECT_FLOAT_EQ(12.3, var_value<T>(new vari_value<T>(12.3)).val())
      << "Failed For T: " << type_name<T>() << std::endl;
  // make sure rvalue var_value can be accepted
  EXPECT_FLOAT_EQ(12.3, var_value<T>(var_value<T>(12.3)).val())
      << "Failed For T: " << type_name<T>() << std::endl;
  // S type is preserved
  EXPECT_FLOAT_EQ(static_cast<S>(3.7), var_value<T>(static_cast<S>(3.7)).val())
      << "Failed For T: " << type_name<T>() << " and S: " << type_name<S>()
      << "\n";
  // Make sure integral types don't hold a nullptr instead of zero.
  EXPECT_FLOAT_EQ(0, var_value<T>(static_cast<S>(0)).val())
      << "Failed For T: " << type_name<T>() << " and S:" << type_name<S>()
      << "\n";
}

template <typename T>
void ctor_overloads_float() {
  ctor_overloads_float_impl<T, double>();
  ctor_overloads_float_impl<T, long double>();
  ctor_overloads_float_impl<T, float>();
  ctor_overloads_float_impl<T, bool>();
  ctor_overloads_float_impl<T, char>();
  ctor_overloads_float_impl<T, int>();
  ctor_overloads_float_impl<T, int16_t>();
  ctor_overloads_float_impl<T, int32_t>();
  ctor_overloads_float_impl<T, unsigned char>();
  ctor_overloads_float_impl<T, unsigned int>();
  ctor_overloads_float_impl<T, uint32_t>();
  ctor_overloads_float_impl<T, uint16_t>();
  ctor_overloads_float_impl<T, size_t>();
  ctor_overloads_float_impl<T, ptrdiff_t>();
}

template <typename EigenMat>
void ctor_overloads_matrix(EigenMat&& xx) {
  using stan::math::var_value;
  using stan::math::vari_value;
  using stan::math::test::type_name;
  using eigen_plain = std::decay_t<stan::plain_type_t<EigenMat>>;

  eigen_plain x = xx;
  // standard constructor
  EXPECT_MATRIX_FLOAT_EQ((x + x).eval(), var_value<eigen_plain>(x + x).val());
  // make sure copy ctor is used rather than casting vari* to unsigned int
  EXPECT_MATRIX_FLOAT_EQ(
      x, var_value<eigen_plain>(new vari_value<eigen_plain>(x)).val());
  // make sure rvalue var_value can be accepted
  EXPECT_MATRIX_FLOAT_EQ(
      x, var_value<eigen_plain>(var_value<eigen_plain>(x)).val());
  // test init_dependent for adj
  auto test_var_x = var_value<eigen_plain>(var_value<eigen_plain>(x));
  test_var_x.vi_->init_dependent();
  EXPECT_MATRIX_FLOAT_EQ(eigen_plain::Ones(x.rows(), x.cols()),
                         test_var_x.adj());
}

template <typename EigenMat>
void ctor_overloads_sparse_matrix(EigenMat&& x) {
  using stan::math::var_value;
  using stan::math::vari_value;
  using stan::math::test::type_name;
  using eigen_plain = std::decay_t<stan::plain_type_t<EigenMat>>;
  using inner_iterator = typename eigen_plain::InnerIterator;
  // standard constructor with eigen expression
  eigen_plain matmul_x = x * x;
  eigen_plain matmul_xx = var_value<eigen_plain>(x * x).val();
  for (int k = 0; k < matmul_x.outerSize(); ++k) {
    for (inner_iterator it(matmul_x, k), iz(matmul_xx, k); it; ++it, ++iz) {
      EXPECT_FLOAT_EQ(iz.value(), it.value());
    }
  }
  const eigen_plain const_matmul_x = x * x;
  eigen_plain const_matmul_xx = var_value<eigen_plain>(const_matmul_x).val();
  for (int k = 0; k < matmul_x.outerSize(); ++k) {
    for (inner_iterator it(const_matmul_x, k), iz(const_matmul_xx, k); it;
         ++it, ++iz) {
      EXPECT_FLOAT_EQ(iz.value(), it.value());
    }
  }

  // make sure rvalue var_value can be accepted
  eigen_plain x_rv = var_value<eigen_plain>(var_value<eigen_plain>(x)).val();
  for (int k = 0; k < x.outerSize(); ++k) {
    for (inner_iterator it(x, k), iz(x_rv, k); it; ++it, ++iz) {
      EXPECT_FLOAT_EQ(iz.value(), it.value());
    }
  }

  // from a vari_value with sparse eigen expression
  eigen_plain x_from_vari
      = var_value<eigen_plain>(new vari_value<eigen_plain>(x * x)).val();
  for (int k = 0; k < matmul_x.outerSize(); ++k) {
    for (inner_iterator it(matmul_x, k), iz(x_from_vari, k); it; ++it, ++iz) {
      EXPECT_FLOAT_EQ(iz.value(), it.value());
    }
  }
  // test inplace addition works
  auto inplace_add_var = var_value<eigen_plain>(new vari_value<eigen_plain>(x));
  eigen_plain test_y = make_sparse_matrix_random(10, 10);
  inplace_add_var.vi_->init_dependent();
  inplace_add_var.adj() += test_y;
  // adjoints sparsity pattern will be pattern of x and test_y for addition
  for (int k = 0; k < x.outerSize(); ++k) {
    for (inner_iterator it(test_y, k), iz(inplace_add_var.adj(), k); iz; ++iz) {
      if (iz.row() == it.row() && iz.col() == it.col()) {
        EXPECT_FLOAT_EQ(iz.value() - 1, it.value());
        ++it;
      } else {
        EXPECT_FLOAT_EQ(iz.value(), 1.0);
      }
    }
  }
}

}  // namespace test
}  // namespace stan
TEST_F(AgradRev, ctorfloatOverloads) {
  stan::test::ctor_overloads_float<float>();
  stan::test::ctor_overloads_float<double>();
  stan::test::ctor_overloads_float<long double>();
}

TEST_F(AgradRev, ctormatrixOverloads) {
  using dense_mat = Eigen::Matrix<double, -1, -1>;
  using sparse_mat = Eigen::SparseMatrix<double>;
  stan::test::ctor_overloads_matrix(dense_mat::Random(10, 10));
  using dense_vec = Eigen::Matrix<double, -1, 1>;
  stan::test::ctor_overloads_matrix(dense_vec::Random(10));
  using dense_row_vec = Eigen::Matrix<double, 1, -1>;
  stan::test::ctor_overloads_matrix(dense_row_vec::Random(10));
  sparse_mat sparse_x = stan::test::make_sparse_matrix_random(10, 10);
  stan::test::ctor_overloads_sparse_matrix(sparse_x);
}

TEST_F(AgradRev, ctorMatrixArenaOverload) {
  using stan::math::arena_matrix;
  using stan::math::var_value;
  arena_matrix<Eigen::MatrixXd> x(Eigen::MatrixXd::Random(5, 5));
  var_value<Eigen::MatrixXd> A(x);
  EXPECT_MATRIX_FLOAT_EQ(A.val(), x);
  const auto& x_ref = x;
  var_value<Eigen::MatrixXd> B(x_ref);
  EXPECT_MATRIX_FLOAT_EQ(B.val(), x);
}

TEST_F(AgradRev, var_matrix_views) {
  using dense_mat = Eigen::Matrix<double, -1, -1>;
  dense_mat A(10, 10);
  for (Eigen::Index i = 0; i < A.size(); ++i) {
    A(i) = i;
  }
  stan::math::var_value<dense_mat> A_v(A);
  auto A_block = A_v.block(1, 1, 3, 3);
  EXPECT_MATRIX_FLOAT_EQ(A_block.val(), A.block(1, 1, 3, 3));
  auto A_transpose = A_v.transpose();
  EXPECT_MATRIX_FLOAT_EQ(A_transpose.val(), A.transpose());
  auto A_row = A_v.row(3);
  EXPECT_MATRIX_FLOAT_EQ(A_row.val(), A.row(3));
  auto A_col = A_v.col(3);
  EXPECT_MATRIX_FLOAT_EQ(A_col.val(), A.col(3));
  auto A_block_row = A_v.block(1, 1, 3, 3).row(1);
  EXPECT_MATRIX_FLOAT_EQ(A_block_row.val(), A.block(1, 1, 3, 3).row(1));
  auto A_rowwise_reverse = A_v.rowwise_reverse();
  EXPECT_MATRIX_FLOAT_EQ(A_rowwise_reverse.val(), A.rowwise().reverse());
  auto A_colwise_reverse = A_v.colwise_reverse();
  EXPECT_MATRIX_FLOAT_EQ(A_colwise_reverse.val(), A.colwise().reverse());
  auto A_rowwise_colwise_reverse = A_v.rowwise_reverse().colwise_reverse();
  EXPECT_MATRIX_FLOAT_EQ(A_rowwise_colwise_reverse.val(),
                         A.rowwise().reverse().colwise().reverse());

  auto A_diagonal = A_v.diagonal();
  EXPECT_MATRIX_FLOAT_EQ(A_diagonal.val(), A.diagonal());

  auto A_coeff1 = A_v(3);
  EXPECT_FLOAT_EQ(A(3), A_coeff1.val());
  auto A_coeff2 = A_v(3, 3);
  EXPECT_FLOAT_EQ(A(3, 3), A_coeff2.val());
  EXPECT_MATRIX_FLOAT_EQ(A, A_v.val());
  for (int i = 0; i < A.size(); ++i) {
    A_v.vi_->adj_(i) = i;
  }
  EXPECT_MATRIX_FLOAT_EQ(A_block.adj(), A_v.adj().block(1, 1, 3, 3));
  EXPECT_MATRIX_FLOAT_EQ(A_transpose.adj(), A_v.adj().transpose());
  EXPECT_MATRIX_FLOAT_EQ(A_row.adj(), A_v.adj().row(3));
  EXPECT_MATRIX_FLOAT_EQ(A_col.adj(), A_v.adj().col(3));
  EXPECT_MATRIX_FLOAT_EQ(A_block_row.adj(), A_v.adj().block(1, 1, 3, 3).row(1));
  EXPECT_MATRIX_FLOAT_EQ(A_rowwise_reverse.adj(),
                         A_v.adj().rowwise().reverse());
  EXPECT_MATRIX_FLOAT_EQ(A_colwise_reverse.adj(),
                         A_v.adj().colwise().reverse());
  EXPECT_MATRIX_FLOAT_EQ(A_rowwise_colwise_reverse.adj(),
                         A_v.adj().rowwise().reverse().colwise().reverse());
  // since new var is made and values propogate back
  A_coeff1.vi_->adj_ = 1;
  A_coeff2.vi_->adj_ = 10;
  auto prev_adj_val1 = A_v.adj()(3);
  auto prev_adj_val2 = A_v.adj()(3, 3);
  stan::math::grad();
  EXPECT_FLOAT_EQ(A_v.adj()(3) - prev_adj_val1, A_coeff1.adj());
  EXPECT_FLOAT_EQ(A_v.adj()(3, 3) - prev_adj_val2, A_coeff2.adj());
}

TEST_F(AgradRev, var_matrix_views_specializations) {
  using dense_mat = Eigen::Matrix<double, -1, -1>;
  dense_mat A(10, 10);
  for (Eigen::Index i = 0; i < A.size(); ++i) {
    A(i) = i;
  }
  stan::math::var_value<dense_mat> A_v(A);
  auto A_toprow = A_v.topRows(1);
  EXPECT_MATRIX_FLOAT_EQ(A_toprow.val(), A.topRows(1));

  auto A_bottomrow = A_v.bottomRows(1);
  EXPECT_MATRIX_FLOAT_EQ(A_bottomrow.val(), A.bottomRows(1));

  auto A_middlerows = A_v.middleRows(3, 2);
  EXPECT_MATRIX_FLOAT_EQ(A_middlerows.val(), A.middleRows(3, 2));

  auto A_leftcols = A_v.leftCols(1);
  EXPECT_MATRIX_FLOAT_EQ(A_leftcols.val(), A.leftCols(1));

  auto A_rightcols = A_v.rightCols(1);
  EXPECT_MATRIX_FLOAT_EQ(A_rightcols.val(), A.rightCols(1));

  auto A_middlecols = A_v.middleCols(3, 2);
  EXPECT_MATRIX_FLOAT_EQ(A_middlecols.val(), A.middleCols(3, 2));

  EXPECT_MATRIX_FLOAT_EQ(A, A_v.val());
  for (int i = 0; i < A.size(); ++i) {
    A_v.vi_->adj_(i) = i;
  }
  EXPECT_MATRIX_FLOAT_EQ(A_toprow.adj(), A_v.adj().topRows(1));
  EXPECT_MATRIX_FLOAT_EQ(A_bottomrow.adj(), A_v.adj().bottomRows(1));
  EXPECT_MATRIX_FLOAT_EQ(A_middlerows.adj(), A_v.adj().middleRows(3, 2));
  EXPECT_MATRIX_FLOAT_EQ(A_leftcols.adj(), A_v.adj().leftCols(1));
  EXPECT_MATRIX_FLOAT_EQ(A_rightcols.adj(), A_v.adj().rightCols(1));
  EXPECT_MATRIX_FLOAT_EQ(A_middlecols.adj(), A_v.adj().middleCols(3, 2));
}

TEST_F(AgradRev, var_matrix_views_const) {
  using dense_mat = Eigen::Matrix<double, -1, -1>;
  dense_mat A(10, 10);
  for (Eigen::Index i = 0; i < A.size(); ++i) {
    A(i) = i;
  }
  stan::math::var_value<dense_mat> A_vv(A);
  const auto& A_v = A_vv;
  auto A_block = A_v.block(1, 1, 3, 3);
  EXPECT_MATRIX_FLOAT_EQ(A_block.val(), A.block(1, 1, 3, 3));
  auto A_transpose = A_v.transpose();
  EXPECT_MATRIX_FLOAT_EQ(A_transpose.val(), A.transpose());
  auto A_row = A_v.row(3);
  EXPECT_MATRIX_FLOAT_EQ(A_row.val(), A.row(3));
  auto A_col = A_v.col(3);
  EXPECT_MATRIX_FLOAT_EQ(A_col.val(), A.col(3));
  auto A_block_row = A_v.block(1, 1, 3, 3).row(1);
  EXPECT_MATRIX_FLOAT_EQ(A_block_row.val(), A.block(1, 1, 3, 3).row(1));
  auto A_rowwise_reverse = A_v.rowwise_reverse();
  EXPECT_MATRIX_FLOAT_EQ(A_rowwise_reverse.val(), A.rowwise().reverse());
  auto A_colwise_reverse = A_v.colwise_reverse();
  EXPECT_MATRIX_FLOAT_EQ(A_colwise_reverse.val(), A.colwise().reverse());
  auto A_rowwise_colwise_reverse = A_v.rowwise_reverse().colwise_reverse();
  EXPECT_MATRIX_FLOAT_EQ(A_rowwise_colwise_reverse.val(),
                         A.rowwise().reverse().colwise().reverse());
  auto A_coeff1 = A_v(3);
  EXPECT_FLOAT_EQ(A(3), A_coeff1.val());
  auto A_coeff2 = A_v(3, 3);
  EXPECT_FLOAT_EQ(A(3, 3), A_coeff2.val());
  EXPECT_MATRIX_FLOAT_EQ(A, A_v.val());
  for (int i = 0; i < A.size(); ++i) {
    A_v.vi_->adj_(i) = i;
  }
  EXPECT_MATRIX_FLOAT_EQ(A_block.adj(), A_v.adj().block(1, 1, 3, 3));
  EXPECT_MATRIX_FLOAT_EQ(A_transpose.adj(), A_v.adj().transpose());
  EXPECT_MATRIX_FLOAT_EQ(A_row.adj(), A_v.adj().row(3));
  EXPECT_MATRIX_FLOAT_EQ(A_col.adj(), A_v.adj().col(3));
  EXPECT_MATRIX_FLOAT_EQ(A_block_row.adj(), A_v.adj().block(1, 1, 3, 3).row(1));
  EXPECT_MATRIX_FLOAT_EQ(A_rowwise_reverse.adj(),
                         A_v.adj().rowwise().reverse());
  EXPECT_MATRIX_FLOAT_EQ(A_colwise_reverse.adj(),
                         A_v.adj().colwise().reverse());
  EXPECT_MATRIX_FLOAT_EQ(A_rowwise_colwise_reverse.adj(),
                         A_v.adj().rowwise().reverse().colwise().reverse());
  // since new var is made and values propogate back
  A_coeff1.vi_->adj_ = 1;
  A_coeff2.vi_->adj_ = 10;
  auto prev_adj_val1 = A_v.adj()(3);
  auto prev_adj_val2 = A_v.adj()(3, 3);
  stan::math::grad();
  EXPECT_FLOAT_EQ(A_v.adj()(3) - prev_adj_val1, A_coeff1.adj());
  EXPECT_FLOAT_EQ(A_v.adj()(3, 3) - prev_adj_val2, A_coeff2.adj());
}

template <typename dense_vec>
void var_vector_views_test() {
  using stan::math::var_value;
  dense_vec A(10);
  for (Eigen::Index i = 0; i < A.size(); ++i) {
    A(i) = i;
  }
  var_value<dense_vec> A_v(A);
  auto A_head = A_v.head(3);
  EXPECT_MATRIX_FLOAT_EQ(A.head(3), A_head.val());
  auto A_transpose = A_v.transpose();
  EXPECT_MATRIX_FLOAT_EQ(A.transpose(), A_transpose.val());
  auto A_tail = A_v.tail(3);
  EXPECT_MATRIX_FLOAT_EQ(A.tail(3), A_tail.val());
  auto A_segment = A_v.segment(3, 5);
  EXPECT_MATRIX_FLOAT_EQ(A.segment(3, 5), A_segment.val());
  auto A_coeff1 = A_v(3);
  EXPECT_FLOAT_EQ(A(3), A_coeff1.val());
  EXPECT_MATRIX_FLOAT_EQ(A, A_v.val());
  for (Eigen::Index i = 0; i < A.size(); ++i) {
    A_v.vi_->adj_(i) = i;
  }
  EXPECT_MATRIX_FLOAT_EQ(A_v.adj().head(3), A_head.adj());
  EXPECT_MATRIX_FLOAT_EQ(A_v.adj().transpose(), A_transpose.adj());
  EXPECT_MATRIX_FLOAT_EQ(A_v.adj().tail(3), A_tail.adj());
  EXPECT_MATRIX_FLOAT_EQ(A_v.adj().segment(3, 5), A_segment.adj());
  // since new var is made and values propogate back
  A_coeff1.vi_->adj_ = 1;
  auto prev_adj_val = A_v.adj()(3);
  stan::math::grad();
  EXPECT_FLOAT_EQ(A_v.adj()(3) - prev_adj_val, A_coeff1.adj());
}

TEST_F(AgradRev, var_vector_views) {
  var_vector_views_test<Eigen::VectorXd>();
  var_vector_views_test<Eigen::RowVectorXd>();
}

template <typename dense_vec>
void var_vector_views_const_test() {
  using stan::math::var_value;
  dense_vec A(10);
  for (Eigen::Index i = 0; i < A.size(); ++i) {
    A(i) = i;
  }
  var_value<dense_vec> A_vv(A);
  const auto& A_v = A_vv;
  auto A_head = A_v.head(3);
  EXPECT_MATRIX_FLOAT_EQ(A.head(3), A_head.val());
  auto A_transpose = A_v.transpose();
  EXPECT_MATRIX_FLOAT_EQ(A.transpose(), A_transpose.val());
  auto A_tail = A_v.tail(3);
  EXPECT_MATRIX_FLOAT_EQ(A.tail(3), A_tail.val());
  auto A_segment = A_v.segment(3, 5);
  EXPECT_MATRIX_FLOAT_EQ(A.segment(3, 5), A_segment.val());
  auto A_coeff1 = A_v(3);
  EXPECT_FLOAT_EQ(A(3), A_coeff1.val());
  EXPECT_MATRIX_FLOAT_EQ(A, A_v.val());
  for (Eigen::Index i = 0; i < A.size(); ++i) {
    A_v.vi_->adj_(i) = i;
  }
  EXPECT_MATRIX_FLOAT_EQ(A_v.adj().head(3), A_head.adj());
  EXPECT_MATRIX_FLOAT_EQ(A_v.adj().transpose(), A_transpose.adj());
  EXPECT_MATRIX_FLOAT_EQ(A_v.adj().tail(3), A_tail.adj());
  EXPECT_MATRIX_FLOAT_EQ(A_v.adj().segment(3, 5), A_segment.adj());
  // since new var is made and values propogate back
  A_coeff1.vi_->adj_ = 1;
  auto prev_adj_val = A_v.adj()(3);
  stan::math::grad();
  EXPECT_FLOAT_EQ(A_v.adj()(3) - prev_adj_val, A_coeff1.adj());
}

TEST_F(AgradRev, var_vector_views_const) {
  var_vector_views_const_test<Eigen::VectorXd>();
  var_vector_views_const_test<Eigen::RowVectorXd>();
}

TEST_F(AgradRev, var_matrix_view_block_from_plain_test) {
  using stan::math::sum;
  using stan::math::var_value;
  Eigen::MatrixXd A(4, 4);
  A << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
  Eigen::MatrixXd B(2, 2);
  B << 2, 3, 4, 5;
  var_value<Eigen::MatrixXd> A_v(A);
  var_value<Eigen::MatrixXd> B_v(B);
  A_v.block(0, 0, 2, 2) = B_v;
  for (Eigen::Index i = 0; i < A_v.size(); ++i) {
    A_v.adj().coeffRef(i) = i;
  }
  stan::math::grad();
  EXPECT_FLOAT_EQ(B_v.adj()(0), 0);
  EXPECT_FLOAT_EQ(B_v.adj()(1), 1);
  EXPECT_FLOAT_EQ(B_v.adj()(2), 4);
  EXPECT_FLOAT_EQ(B_v.adj()(3), 5);
}

/**
 * Tests that views of a var<Matrix> receive the adjoints of the original
 * matrix.
 */
TEST_F(AgradRev, var_matrix_view) {
  using stan::math::sum;
  using stan::math::var_value;
  Eigen::MatrixXd A(4, 4);
  A << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
  var_value<Eigen::MatrixXd> A_v(A);
  auto A_v_block = A_v.block(1, 1, 3, 3);
  auto A_v_transpose = A_v.transpose();
  auto A_v_row = A_v.row(3);
  auto A_v_col = A_v.col(3);
  auto A_v_block_row = A_v.block(1, 1, 3, 3).row(1);
  auto A_v_rowwise_reverse = A_v.rowwise_reverse();
  auto A_v_colwise_reverse = A_v.colwise_reverse();
  auto A_v_rowwise_colwise_reverse = A_v.rowwise_reverse().colwise_reverse();
  // NOTE: Coefficient references make a new var.
  auto A_v_coeff1 = A_v.coeff(5);
  auto A_v_coeff2 = A_v.coeff(1, 2);
  A_v.block(0, 0, 3, 3) = A_v.block(1, 1, 3, 3);
  stan::math::sum(A_v).grad();
  Eigen::MatrixXd deriv(4, 4);
  deriv << 0, 0, 0, 1, 0, 1, 1, 2, 0, 1, 1, 2, 1, 2, 2, 2;

  EXPECT_MATRIX_FLOAT_EQ(A_v.adj(), deriv);
  EXPECT_MATRIX_FLOAT_EQ(A_v_block.val(), A_v.val().block(1, 1, 3, 3));
  EXPECT_MATRIX_FLOAT_EQ(A_v_block.adj(), A_v.adj().block(1, 1, 3, 3));

  EXPECT_MATRIX_FLOAT_EQ(A_v_transpose.val(), A_v.val().transpose());
  EXPECT_MATRIX_FLOAT_EQ(A_v_transpose.adj(), A_v.adj().transpose());

  EXPECT_MATRIX_FLOAT_EQ(A_v_row.val(), A_v.val().row(3));
  EXPECT_MATRIX_FLOAT_EQ(A_v_row.adj(), A_v.adj().row(3));

  EXPECT_MATRIX_FLOAT_EQ(A_v_col.val(), A_v.val().col(3));
  EXPECT_MATRIX_FLOAT_EQ(A_v_col.adj(), A_v.adj().col(3));

  EXPECT_MATRIX_FLOAT_EQ(A_v_block_row.val(),
                         A_v.val().block(1, 1, 3, 3).row(1));
  EXPECT_MATRIX_FLOAT_EQ(A_v_block_row.adj(),
                         A_v.adj().block(1, 1, 3, 3).row(1));

  EXPECT_MATRIX_FLOAT_EQ(A_v_rowwise_reverse.val(),
                         A_v.val().rowwise().reverse());
  EXPECT_MATRIX_FLOAT_EQ(A_v_rowwise_reverse.adj(),
                         A_v.adj().rowwise().reverse());

  EXPECT_MATRIX_FLOAT_EQ(A_v_colwise_reverse.val(),
                         A_v.val().colwise().reverse());
  EXPECT_MATRIX_FLOAT_EQ(A_v_colwise_reverse.adj(),
                         A_v.adj().colwise().reverse());

  EXPECT_MATRIX_FLOAT_EQ(A_v_rowwise_colwise_reverse.val(),
                         A_v.val().rowwise().reverse().colwise().reverse());
  EXPECT_MATRIX_FLOAT_EQ(A_v_rowwise_colwise_reverse.adj(),
                         A_v.adj().rowwise().reverse().colwise().reverse());

  EXPECT_FLOAT_EQ(A_v_coeff1.val(), A_v.val().coeff(5));
  EXPECT_FLOAT_EQ(A_v_coeff1.adj(), 0);

  EXPECT_FLOAT_EQ(A_v_coeff2.val(), A_v.val().coeff(1, 2));
  EXPECT_FLOAT_EQ(A_v_coeff2.adj(), 0);
}

TEST_F(AgradRev, var_matrix_view_const) {
  using stan::math::sum;
  using stan::math::var_value;
  Eigen::MatrixXd A(4, 4);
  A << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
  var_value<Eigen::MatrixXd> A_vv(A);
  const auto& A_v = A_vv;
  auto A_v_block = A_v.block(1, 1, 3, 3);
  auto A_v_transpose = A_v.transpose();
  auto A_v_row = A_v.row(3);
  auto A_v_col = A_v.col(3);
  auto A_v_block_row = A_v.block(1, 1, 3, 3).row(1);
  auto A_v_rowwise_reverse = A_v.rowwise_reverse();
  auto A_v_colwise_reverse = A_v.colwise_reverse();
  auto A_v_rowwise_colwise_reverse = A_v.rowwise_reverse().colwise_reverse();
  // NOTE: Coefficient references make a new var.
  auto A_v_coeff1 = A_v.coeff(5);
  auto A_v_coeff2 = A_v.coeff(1, 2);
  A_v.block(0, 0, 3, 3) = A_v.block(1, 1, 3, 3);
  stan::math::sum(A_v).grad();
  Eigen::MatrixXd deriv(4, 4);
  deriv << 0, 0, 0, 1, 0, 1, 1, 2, 0, 1, 1, 2, 1, 2, 2, 2;

  EXPECT_MATRIX_FLOAT_EQ(A_v.adj(), deriv);
  EXPECT_MATRIX_FLOAT_EQ(A_v_block.val(), A_v.val().block(1, 1, 3, 3));
  EXPECT_MATRIX_FLOAT_EQ(A_v_block.adj(), A_v.adj().block(1, 1, 3, 3));

  EXPECT_MATRIX_FLOAT_EQ(A_v_transpose.val(), A_v.val().transpose());
  EXPECT_MATRIX_FLOAT_EQ(A_v_transpose.adj(), A_v.adj().transpose());

  EXPECT_MATRIX_FLOAT_EQ(A_v_row.val(), A_v.val().row(3));
  EXPECT_MATRIX_FLOAT_EQ(A_v_row.adj(), A_v.adj().row(3));

  EXPECT_MATRIX_FLOAT_EQ(A_v_col.val(), A_v.val().col(3));
  EXPECT_MATRIX_FLOAT_EQ(A_v_col.adj(), A_v.adj().col(3));

  EXPECT_MATRIX_FLOAT_EQ(A_v_block_row.val(),
                         A_v.val().block(1, 1, 3, 3).row(1));
  EXPECT_MATRIX_FLOAT_EQ(A_v_block_row.adj(),
                         A_v.adj().block(1, 1, 3, 3).row(1));

  EXPECT_MATRIX_FLOAT_EQ(A_v_rowwise_reverse.val(),
                         A_v.val().rowwise().reverse());
  EXPECT_MATRIX_FLOAT_EQ(A_v_rowwise_reverse.adj(),
                         A_v.adj().rowwise().reverse());

  EXPECT_MATRIX_FLOAT_EQ(A_v_colwise_reverse.val(),
                         A_v.val().colwise().reverse());
  EXPECT_MATRIX_FLOAT_EQ(A_v_colwise_reverse.adj(),
                         A_v.adj().colwise().reverse());

  EXPECT_MATRIX_FLOAT_EQ(A_v_rowwise_colwise_reverse.val(),
                         A_v.val().rowwise().reverse().colwise().reverse());
  EXPECT_MATRIX_FLOAT_EQ(A_v_rowwise_colwise_reverse.adj(),
                         A_v.adj().rowwise().reverse().colwise().reverse());

  EXPECT_FLOAT_EQ(A_v_coeff1.val(), A_v.val().coeff(5));
  EXPECT_FLOAT_EQ(A_v_coeff1.adj(), 0);

  EXPECT_FLOAT_EQ(A_v_coeff2.val(), A_v.val().coeff(1, 2));
  EXPECT_FLOAT_EQ(A_v_coeff2.adj(), 0);
}

TEST_F(AgradRev, var_matrix_view_assignment) {
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  Eigen::MatrixXd A(4, 4);
  A << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
  var_value<Eigen::MatrixXd> A_v(A);
  var_value<Eigen::MatrixXd> A_v_block = A_v.block(1, 1, 3, 3);
  var_value<Eigen::MatrixXd> A_v_transpose = A_v.transpose();
  var_value<Eigen::RowVectorXd> A_v_row = A_v.row(3);
  var_value<Eigen::VectorXd> A_v_col = A_v.col(3);
  var_value<Eigen::RowVectorXd> A_v_block_row = A_v.block(1, 1, 3, 3).row(1);
  var_value<Eigen::MatrixXd> A_v_rowwise_reverse = A_v.rowwise_reverse();
  var_value<Eigen::MatrixXd> A_v_colwise_reverse = A_v.colwise_reverse();
  var_value<Eigen::MatrixXd> A_v_rowwise_colwise_reverse
      = A_v.rowwise_reverse().colwise_reverse();
  var A_v_coeff1 = A_v.coeff(5);
  var A_v_coeff2 = A_v.coeff(1, 2);
  A_v.block(0, 0, 3, 3) = A_v.block(1, 1, 3, 3);
  // Checks adjoints from all assigned slices are propogated upwards
  var b_v = stan::math::sum(A_v_block) + stan::math::sum(A_v_transpose)
            + stan::math::sum(A_v_row) + stan::math::sum(A_v_col)
            + stan::math::sum(A_v_block_row)
            + stan::math::sum(A_v_rowwise_reverse)
            + stan::math::sum(A_v_colwise_reverse)
            + stan::math::sum(A_v_rowwise_colwise_reverse);
  b_v.grad();
  Eigen::MatrixXd deriv(4, 4);
  deriv << 4, 4, 4, 5, 4, 5, 5, 6, 4, 6, 6, 7, 5, 6, 6, 7;
  EXPECT_MATRIX_FLOAT_EQ(A_v.adj(), deriv);
}

TEST_F(AgradRev, var_matrix_view_assignment_const) {
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  Eigen::MatrixXd A(4, 4);
  A << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
  var_value<Eigen::MatrixXd> A_vv(A);
  const auto& A_v = A_vv;
  var_value<Eigen::MatrixXd> A_v_block = A_v.block(1, 1, 3, 3);
  var_value<Eigen::MatrixXd> A_v_transpose = A_v.transpose();
  var_value<Eigen::RowVectorXd> A_v_row = A_v.row(3);
  var_value<Eigen::VectorXd> A_v_col = A_v.col(3);
  var_value<Eigen::RowVectorXd> A_v_block_row = A_v.block(1, 1, 3, 3).row(1);
  var_value<Eigen::MatrixXd> A_v_rowwise_reverse = A_v.rowwise_reverse();
  var_value<Eigen::MatrixXd> A_v_colwise_reverse = A_v.colwise_reverse();
  var_value<Eigen::MatrixXd> A_v_rowwise_colwise_reverse
      = A_v.rowwise_reverse().colwise_reverse();
  var A_v_coeff1 = A_v.coeff(5);
  var A_v_coeff2 = A_v.coeff(1, 2);
  A_v.block(0, 0, 3, 3) = A_v.block(1, 1, 3, 3);
  // Checks adjoints from all assigned slices are propogated upwards
  var b_v = stan::math::sum(A_v_block) + stan::math::sum(A_v_transpose)
            + stan::math::sum(A_v_row) + stan::math::sum(A_v_col)
            + stan::math::sum(A_v_block_row)
            + stan::math::sum(A_v_rowwise_reverse)
            + stan::math::sum(A_v_colwise_reverse)
            + stan::math::sum(A_v_rowwise_colwise_reverse);
  b_v.grad();
  Eigen::MatrixXd deriv(4, 4);
  deriv << 4, 4, 4, 5, 4, 5, 5, 6, 4, 6, 6, 7, 5, 6, 6, 7;
  EXPECT_MATRIX_FLOAT_EQ(A_v.adj(), deriv);
}

TEST_F(AgradRev, var_matrix_view_eval) {
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  Eigen::MatrixXd A(4, 4);
  A << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
  var_value<Eigen::MatrixXd> A_v(A);
  auto A_v_block = A_v.block(1, 1, 3, 3).eval();
  auto A_v_transpose = A_v.transpose().eval();
  auto A_v_row = A_v.row(3).eval();
  auto A_v_col = A_v.col(3).eval();
  auto A_v_block_row = A_v.block(1, 1, 3, 3).row(1).eval();
  auto A_v_rowwise_reverse = A_v.rowwise_reverse().eval();
  auto A_v_colwise_reverse = A_v.colwise_reverse().eval();
  auto A_v_rowwise_colwise_reverse
      = A_v.rowwise_reverse().colwise_reverse().eval();
  // NOTE: Coefficient references make a new var.
  auto A_v_coeff1 = A_v.coeff(5);
  auto A_v_coeff2 = A_v.coeff(1, 2);
  A_v.block(0, 0, 3, 3) = A_v.block(1, 1, 3, 3);
  // Checks adjoints from all assigned slices are propogated upwards
  var b_v = stan::math::sum(A_v_block) + stan::math::sum(A_v_transpose)
            + stan::math::sum(A_v_row) + stan::math::sum(A_v_col)
            + stan::math::sum(A_v_block_row)
            + stan::math::sum(A_v_rowwise_reverse)
            + stan::math::sum(A_v_colwise_reverse)
            + stan::math::sum(A_v_rowwise_colwise_reverse);
  b_v.grad();
  Eigen::MatrixXd deriv(4, 4);
  deriv << 4, 4, 4, 5, 4, 5, 5, 6, 4, 6, 6, 7, 5, 6, 6, 7;
  EXPECT_MATRIX_FLOAT_EQ(A_v.adj(), deriv);
}

TEST_F(AgradRev, var_matrix_view_block_plain_assignment) {
  Eigen::MatrixXd A(4, 4);
  A << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
  stan::math::var_value<Eigen::MatrixXd> A_v(A);
  stan::math::var_value<Eigen::MatrixXd> B_v = A_v.block(1, 1, 3, 3);
  stan::math::sum(B_v).grad();
  Eigen::MatrixXd deriv(4, 4);
  deriv << 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1;
  EXPECT_MATRIX_FLOAT_EQ(A_v.adj(), deriv);
}

TEST_F(AgradRev, var_matrix_view_transpose_plain_assignment) {
  Eigen::MatrixXd A(4, 4);
  A << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
  stan::math::var_value<Eigen::MatrixXd> A_v(A);
  stan::math::var_value<Eigen::MatrixXd> B_v = A_v.transpose();
  stan::math::sum(B_v).grad();
  Eigen::MatrixXd deriv(4, 4);
  deriv << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;
  EXPECT_MATRIX_FLOAT_EQ(A_v.adj(), deriv);
}

TEST_F(AgradRev, var_matrix_view_row_plain_assignment) {
  Eigen::MatrixXd A(4, 4);
  A << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
  stan::math::var_value<Eigen::MatrixXd> A_v(A);
  stan::math::var_value<Eigen::RowVectorXd> B_v = A_v.row(3);
  stan::math::var b_v = A_v(1) + A_v(1, 1) + stan::math::sum(B_v);
  b_v.grad();
  Eigen::MatrixXd deriv(4, 4);
  deriv << 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1;
  EXPECT_MATRIX_FLOAT_EQ(A_v.adj(), deriv);
}

TEST_F(AgradRev, var_matrix_array) {
  Eigen::MatrixXd A(4, 4);
  A << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
  stan::math::var_value<Eigen::MatrixXd> A_v(A);
  Eigen::Array<double, -1, -1> B_v = A_v.array().val();
}

TEST_F(AgradRev, var_matrix_vector_rowvector_assign) {
  Eigen::RowVectorXd A = Eigen::RowVectorXd(4);
  stan::math::var_value<Eigen::VectorXd> A_v(A);
}

TEST_F(AgradRev, a_eq_x) {
  stan::math::var a = 5.0;
  EXPECT_FLOAT_EQ(5.0, a.val());
}

TEST_F(AgradRev, a_of_x) {
  stan::math::var a(6.0);
  EXPECT_FLOAT_EQ(6.0, a.val());
}

TEST_F(AgradRev, a__a_eq_x) {
  stan::math::var a;
  a = 7.0;
  EXPECT_FLOAT_EQ(7.0, a.val());
}

TEST_F(AgradRev, eq_a) {
  stan::math::var a = 5.0;
  stan::math::var f = a;
  std::vector<stan::math::var> x{a};
  std::vector<double> dx;
  f.grad(x, dx);
  EXPECT_FLOAT_EQ(1.0, dx[0]);
}

TEST_F(AgradRev, a_ostream) {
  stan::math::var a = 6.0;
  std::ostringstream os;

  os << a;
  EXPECT_EQ("6", os.str());

  os.str("");
  a = 10.5;
  os << a;
  EXPECT_EQ("10.5", os.str());
}

TEST_F(AgradRev, smart_ptrs) {
  stan::math::var a = 2.0;
  EXPECT_FLOAT_EQ(2.0, (*a).val_);
  EXPECT_FLOAT_EQ(2.0, a->val_);

  EXPECT_FLOAT_EQ(2.0, (*a.vi_).val_);
  EXPECT_FLOAT_EQ(2.0, a.vi_->val_);
}

TEST_F(AgradRev, stackAllocation) {
  using stan::math::var;
  using stan::math::vari;

  vari ai(1.0);
  vari bi(2.0);

  var a(&ai);
  var b(&bi);

  std::vector<stan::math::var> x{a, b};
  var f = a * b;

  std::vector<double> g;
  f.grad(x, g);

  EXPECT_EQ(2U, g.size());
  EXPECT_FLOAT_EQ(2.0, g[0]);
  EXPECT_FLOAT_EQ(1.0, g[1]);
}

TEST_F(AgradRev, print) {
  using stan::math::var;

  std::ostringstream output;
  std::string str;

  var initialized_var(0);
  output << initialized_var;
  str = output.str();
  EXPECT_STREQ("0", output.str().c_str());

  output.clear();
  output.str("");
  var uninitialized_var;
  output << uninitialized_var;
  str = output.str();
  EXPECT_STREQ("uninitialized", output.str().c_str());
}

// should really be doing this test with a mock object using ctor
// vari_(double, bool);  as in:
//
// struct nostack_test_vari : public stan::math::vari {
//   nostack_test_vari(double x)
//   : stan::math::vari(x, false) {
//   }
//   void chain() {
//     // no op on the chain
//   }
// };

// struct both_test_vari : public stan::math::vari {
//   both_test_vari(vari* vi, vari* bi) {

//   }
// };

// var foo(var y, var z) {
//   return y *
// }

TEST_F(AgradRev, basicGradient1) {
  using stan::math::recover_memory;

  for (int i = 0; i < 100; ++i) {
    gradable g = setup_simple();
    g.test();
    recover_memory();
  }
}

TEST_F(AgradRev, basicGradient2) {
  using stan::math::recover_memory;

  for (int i = 0; i < 100; ++i) {
    gradable g = setup_quad_form();
    g.test();
    recover_memory();
  }
}

TEST_F(AgradRev, nestedGradient1) {
  using stan::math::recover_memory;
  using stan::math::recover_memory_nested;
  using stan::math::start_nested;

  gradable g0 = setup_simple();

  start_nested();
  gradable g1 = setup_quad_form();
  g1.test();
  recover_memory_nested();

  start_nested();
  gradable g2 = setup_simple();
  g2.test();
  recover_memory_nested();

  g0.test();
  recover_memory();
}

TEST_F(AgradRev, nestedGradient2) {
  using stan::math::recover_memory;
  using stan::math::recover_memory_nested;
  using stan::math::start_nested;

  gradable g0 = setup_quad_form();

  start_nested();
  gradable g1 = setup_simple();
  g1.test();
  recover_memory_nested();

  start_nested();
  gradable g2 = setup_quad_form();
  g2.test();
  recover_memory_nested();

  g0.test();
  recover_memory();
}

TEST_F(AgradRev, nestedGradient3) {
  using stan::math::recover_memory;
  using stan::math::recover_memory_nested;
  using stan::math::start_nested;

  start_nested();
  gradable g1 = setup_simple();
  start_nested();
  gradable g2 = setup_quad_form();
  start_nested();
  gradable g3 = setup_quad_form();
  start_nested();
  gradable g4 = setup_simple();
  g4.test();
  recover_memory_nested();
  g3.test();
  recover_memory_nested();
  g2.test();
  recover_memory_nested();
  g1.test();
  recover_memory_nested();
  recover_memory();
}

TEST_F(AgradRev, grad) {
  stan::math::var a = 5.0;
  stan::math::var b = 10.0;
  stan::math::var f = a * b + a;

  EXPECT_NO_THROW(f.grad()) << "testing the grad function with no args";

  EXPECT_FLOAT_EQ(5.0, a.val());
  EXPECT_FLOAT_EQ(10.0, b.val());
  EXPECT_FLOAT_EQ(55.0, f.val());

  EXPECT_FLOAT_EQ(1.0, f.adj());
  EXPECT_FLOAT_EQ(11.0, a.adj());
  EXPECT_FLOAT_EQ(5.0, b.adj());
}

TEST_F(AgradRev, matrix_compile_time_conversions) {
  using stan::math::var_value;
  Eigen::VectorXd colvec_vals = Eigen::VectorXd::Random(5);
  Eigen::RowVectorXd rowvec_vals = Eigen::VectorXd::Random(5);
  Eigen::Matrix<double, 1, 1> x11_vals
      = Eigen::Matrix<double, 1, 1>::Random(1, 1);
  var_value<Eigen::Matrix<double, -1, 1>> colvec = colvec_vals;
  var_value<Eigen::Matrix<double, 1, -1>> rowvec = rowvec_vals;
  colvec = rowvec;
  EXPECT_MATRIX_FLOAT_EQ(colvec.val().transpose(), rowvec.val());
  var_value<Eigen::Matrix<double, 1, 1>> x11 = x11_vals;
  colvec = x11;
  rowvec = x11;
  EXPECT_MATRIX_FLOAT_EQ(colvec.val(), rowvec.val());
  EXPECT_MATRIX_FLOAT_EQ(x11.val(), rowvec.val());
}
