#include <stan/math/mix.hpp>
#include <gtest/gtest.h>

using Eigen::Dynamic;
using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;

using stan::math::hessian;
using stan::math::hessian_block_diag;

namespace hessian_block_diag_test {

// fun1(x, y) = x^2 * y + 3 * y^2 (block-diagonal with full 2x2 block)
struct fun1 {
  template <typename T>
  inline T operator()(const Matrix<T, Dynamic, 1>& x) const {
    return x(0) * x(0) * x(1) + 3.0 * x(1) * x(1);
  }
};

// fun2(x, y) = x^2 * y + 3*y^2 + 5*x*y + sin(x) (full Hessian)
struct fun2 {
  template <typename T>
  inline T operator()(const Matrix<T, Dynamic, 1>& x) const {
    using std::sin;
    return x(0) * x(0) * x(1) + 3.0 * x(1) * x(1) + 5.0 * x(0) * x(1)
           + sin(x(0));
  }
};

// exp_diag(x, y) = exp(2x) + exp(y) (diagonal Hessian)
struct exp_diag {
  template <typename T>
  inline T operator()(const Matrix<T, Dynamic, 1>& x) const {
    using stan::math::exp;
    return exp(2 * x(0)) + exp(x(1));
  }
};

// one_arg(x) = x^3 (1x1 Hessian)
struct one_arg {
  template <typename T>
  inline T operator()(const Matrix<T, Dynamic, 1>& x) const {
    return stan::math::pow(x(0), 3);
  }
};

// block_fun: two independent 2x2 blocks on a 4-dimensional input
struct block_fun {
  template <typename T>
  inline T operator()(const Matrix<T, Dynamic, 1>& x) const {
    using stan::math::exp;
    using std::sin;
    // block1: x0^2 + sin(x1)
    // block2: exp(x2) * x3
    return x(0) * x(0) + sin(x(1)) + exp(x(2)) * x(3);
  }
};

TEST(MixFunctor, HessianBlockDiagFun1FullBlock) {
  VectorXd x{{2.0, -3.0}};
  // compute block-diagonal Hessian with block size = full dimension
  Eigen::SparseMatrix<double> H_block = hessian_block_diag(fun1{}, x, 2);

  // compute full dense Hessian
  double fx;
  VectorXd grad;
  MatrixXd H_full;
  hessian(fun1{}, x, fx, grad, H_full);
  EXPECT_EQ(H_full.rows(), H_block.rows());
  EXPECT_EQ(H_full.cols(), H_block.cols());
  for (int i = 0; i < H_full.rows(); ++i) {
    for (int j = 0; j < H_full.cols(); ++j) {
      EXPECT_NEAR(H_full(i, j), H_block.coeff(i, j), 1e-12);
    }
  }
}

TEST(MixFunctor, HessianBlockDiagFun2FullBlock) {
  VectorXd x{{13.0, -4.0}};
  Eigen::SparseMatrix<double> H_block = hessian_block_diag(fun2{}, x, 2);

  double fx;
  VectorXd grad;
  MatrixXd H_full;
  hessian(fun2{}, x, fx, grad, H_full);

  MatrixXd H_block_dense = MatrixXd(H_block);
  EXPECT_EQ(2, H_block_dense.rows());
  EXPECT_EQ(2, H_block_dense.cols());
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      EXPECT_NEAR(H_full(i, j), H_block_dense(i, j), 1e-12);
    }
  }
}

TEST(MixFunctor, HessianBlockDiagExpDiagBlock1) {
  VectorXd x{{2.0, -1.0}};
  Eigen::SparseMatrix<double> H_block = hessian_block_diag(exp_diag{}, x, 1);

  double fx;
  VectorXd grad;
  MatrixXd H_full;
  hessian(exp_diag{}, x, fx, grad, H_full);

  EXPECT_EQ(H_full.rows(), H_block.rows());
  EXPECT_EQ(H_full.cols(), H_block.cols());
  // Only diagonal entries should be captured
  for (int i = 0; i < H_full.rows(); ++i) {
    for (int j = 0; j < H_full.cols(); ++j) {
      if (i == j) {
        EXPECT_NEAR(H_full(i, j), H_block.coeff(i, j), 1e-12);
      } else {
        EXPECT_FLOAT_EQ(0.0, H_block.coeff(i, j));
      }
    }
  }
}

TEST(MixFunctor, HessianBlockDiagOneArg) {
  VectorXd x{{8.0}};

  Eigen::SparseMatrix<double> H_block = hessian_block_diag(one_arg{}, x, 1);

  double fx;
  VectorXd grad;
  MatrixXd H_full;
  hessian(one_arg{}, x, fx, grad, H_full);

  EXPECT_EQ(1, H_block.rows());
  EXPECT_EQ(1, H_block.cols());
  EXPECT_NEAR(H_full(0, 0), H_block.coeff(0, 0), 1e-12);
}

TEST(MixFunctor, HessianBlockDiagBlockFunMultiBlock) {
  VectorXd x{{1.5, -0.5, 0.7, 2.3}};
  // two blocks of size 2
  Eigen::SparseMatrix<double> H_block = hessian_block_diag(block_fun{}, x, 2);
  double fx;
  VectorXd grad;
  MatrixXd H_full;
  hessian(block_fun{}, x, fx, grad, H_full);

  MatrixXd H_block_dense = MatrixXd(H_block);
  // block 0: rows/cols [0,1]
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      EXPECT_NEAR(H_full(i, j), H_block_dense(i, j), 1e-12);
    }
  }
  // block 1: rows/cols [2,3]
  for (int i = 2; i < 4; ++i) {
    for (int j = 2; j < 4; ++j) {
      EXPECT_NEAR(H_full(i, j), H_block_dense(i, j), 1e-12);
    }
  }
  // off-block entries should be zero
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      if (!((i < 2 && j < 2) || (i >= 2 && j >= 2))) {
        EXPECT_FLOAT_EQ(0.0, H_block_dense(i, j));
      }
    }
  }
}

}  // namespace hessian_block_diag_test
