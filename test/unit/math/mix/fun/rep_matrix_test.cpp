#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, repMatrix) {
  // y is scalar
  auto f = [](int m, int n) {
    return [=](const auto& y) { return stan::math::rep_matrix(y, m, n); };
  };

  // y is row vector or column vector
  auto g = [](int k) {
    return [=](const auto& y) { return stan::math::rep_matrix(y, k); };
  };

  double y = 3;
  stan::test::expect_ad(f(0, 0), y);
  stan::test::expect_ad(f(1, 1), y);
  stan::test::expect_ad(f(2, 3), y);

  // illegal arguments---test throw
  stan::test::expect_ad(f(-2, -1), y);

  Eigen::VectorXd a(3);
  a << 3, 3, 3;
  stan::test::expect_ad(g(2), a);

  Eigen::RowVectorXd b(2);
  b << 2, 2;
  stan::test::expect_ad(g(3), b);
}

TEST(MathMixMatFun, repVarMatrix) {
  using stan::math::rep_matrix;
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  auto x_var = var(1.0);
  auto x
      = rep_matrix<var_value<Eigen::Matrix<double, -1, -1>>>(x_var, 5, 5);
  auto x_sum = sum(x);
  x_sum.grad();

  EXPECT_EQ(x_sum.val(), 25.0);
  EXPECT_EQ(x_sum.adj(), 1.0);
  EXPECT_MATRIX_EQ(x.val(), Eigen::MatrixXd::Ones(5, 5));
  EXPECT_MATRIX_EQ(x.adj(), Eigen::MatrixXd::Ones(5, 5));
  EXPECT_EQ(x_var.val(), 1.0);
  EXPECT_EQ(x_var.adj(), 25.0);
}
TEST(MathMixMatFun, repVarMatrixArithmetic) {
  using stan::math::rep_matrix;
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  double x_dbl = 1.0;
  auto x
      = rep_matrix<var_value<Eigen::Matrix<double, -1, -1>>>(x_dbl, 5, 5);
  auto x_sum = sum(x);
  x_sum.grad();

  EXPECT_EQ(x_sum.val(), 25.0);
  EXPECT_EQ(x_sum.adj(), 1.0);
  EXPECT_MATRIX_EQ(x.val(), Eigen::MatrixXd::Ones(5, 5));
  EXPECT_MATRIX_EQ(x.adj(), Eigen::MatrixXd::Ones(5, 5));
}

TEST(MathMixMatFun, repVarMatrixVec) {
  using stan::math::rep_matrix;
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  var_value<Eigen::VectorXd> x_var(Eigen::VectorXd::Ones(5));
  auto x = rep_matrix<var_value<Eigen::Matrix<double, -1, -1>>>(x_var, 5);
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < x.rows(); ++i) {
      x.adj()(i, j) = i;
    }
  }
  auto x_sum = sum(x);
  x_sum.grad();
  EXPECT_EQ(x_sum.val(), 25.0);
  EXPECT_EQ(x_sum.adj(), 1.0);
  Eigen::VectorXd expected_cols(5);
  expected_cols << 1, 2, 3, 4, 5;
  auto expected_adjs = expected_cols.replicate(1, 5).eval();
  EXPECT_MATRIX_EQ(x.val(), Eigen::MatrixXd::Ones(5, 5));
  EXPECT_MATRIX_EQ(x.adj(), expected_adjs);
  EXPECT_MATRIX_EQ(x_var.val(), Eigen::VectorXd::Ones(5));
  Eigen::VectorXd expected_x_var_adjs(5);
  expected_x_var_adjs << 15, 15, 15, 15, 15;
  EXPECT_MATRIX_EQ(x_var.adj(), expected_x_var_adjs);
}

TEST(MathMixMatFun, repVarMatrixVecArithmetic) {
  using stan::math::rep_matrix;
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  auto x_dbl  = Eigen::VectorXd::Ones(5).eval();
  auto x = rep_matrix<var_value<Eigen::Matrix<double, -1, -1>>>(x_dbl, 5);
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < x.rows(); ++i) {
      x.adj()(i, j) = i;
    }
  }
  auto x_sum = sum(x);
  x_sum.grad();
  EXPECT_EQ(x_sum.val(), 25.0);
  EXPECT_EQ(x_sum.adj(), 1.0);
  Eigen::VectorXd expected_cols(5);
  expected_cols << 1, 2, 3, 4, 5;
  auto expected_adjs = expected_cols.replicate(1, 5).eval();
  EXPECT_MATRIX_EQ(x.val(), Eigen::MatrixXd::Ones(5, 5));
  EXPECT_MATRIX_EQ(x.adj(), expected_adjs);
}

TEST(MathMixMatFun, repVarMatrixRowVec) {
  using stan::math::rep_matrix;
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  var_value<Eigen::RowVectorXd> x_var(Eigen::RowVectorXd::Ones(5));
  auto x = rep_matrix<var_value<Eigen::Matrix<double, -1, -1>>>(x_var, 5);
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < x.rows(); ++i) {
      x.adj()(i, j) = i;
    }
  }
  auto x_sum = sum(x);
  x_sum.grad();

  EXPECT_EQ(x_sum.val(), 25.0);
  EXPECT_EQ(x_sum.adj(), 1.0);
  Eigen::VectorXd expected_cols(5);
  expected_cols << 1, 2, 3, 4, 5;
  auto expected_adjs = expected_cols.replicate(1, 5).eval();
  EXPECT_MATRIX_EQ(x.val(), Eigen::MatrixXd::Ones(5, 5));
  EXPECT_MATRIX_EQ(x.adj(), expected_adjs);
  EXPECT_MATRIX_EQ(x_var.val(), Eigen::RowVectorXd::Ones(5));
  Eigen::RowVectorXd expected_x_var_adjs(5);
  expected_x_var_adjs << 5, 10, 15, 20, 25;
  EXPECT_MATRIX_EQ(x_var.adj(), expected_x_var_adjs);
}

TEST(MathMixMatFun, repVarMatrixRowVecArithmetic) {
  using stan::math::rep_matrix;
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  auto x_dbl = Eigen::RowVectorXd::Ones(5).eval();
  auto x = rep_matrix<var_value<Eigen::Matrix<double, -1, -1>>>(x_dbl, 5);
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < x.rows(); ++i) {
      x.adj()(i, j) = i;
    }
  }
  auto x_sum = sum(x);
  x_sum.grad();

  EXPECT_EQ(x_sum.val(), 25.0);
  EXPECT_EQ(x_sum.adj(), 1.0);
  Eigen::VectorXd expected_cols(5);
  expected_cols << 1, 2, 3, 4, 5;
  auto expected_adjs = expected_cols.replicate(1, 5).eval();
  EXPECT_MATRIX_EQ(x.val(), Eigen::MatrixXd::Ones(5, 5));
  EXPECT_MATRIX_EQ(x.adj(), expected_adjs);
}
