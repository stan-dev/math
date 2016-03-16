#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>

struct chol_functor {
  int i, j, K;
  chol_functor(int i_, int j_, int K_) : i(i_), j(j_), K(K_) { }
  template <typename T>
  T operator()(Eigen::Matrix<T, -1, 1> x) const {
    using stan::math::cov_matrix_constrain;
    using stan::math::cholesky_decompose;
    T lp(0.0);
    Eigen::Matrix<T, -1, -1> x_c = cov_matrix_constrain(x, K, lp);
    Eigen::Matrix<T, -1, -1> L = cholesky_decompose(x_c);
    lp += L(i, j);
    return lp;
  }
};

void test_hessians(int size) {
  std::vector<std::vector<chol_functor> > functowns;
  std::vector<std::vector<Eigen::Matrix<double, -1, 1> > > grads_ad;
  std::vector<std::vector<Eigen::Matrix<double, -1, 1> > > grads_fd;
  std::vector<std::vector<Eigen::Matrix<double, -1, -1> > > hess_ad;
  std::vector<std::vector<Eigen::Matrix<double, -1, -1> > > hess_fd;
  Eigen::Matrix<double, -1, -1> evals_ad(size,size);
  Eigen::Matrix<double, -1, -1> evals_fd(size,size);
  functowns.resize(size);
  grads_ad.resize(size);
  grads_fd.resize(size);
  hess_ad.resize(size);
  hess_fd.resize(size);

  for (int i = 0; i < size; ++i)
    for (int j = 0; j < size; ++j) {
      functowns[i].push_back(chol_functor(i, j, size));
      grads_fd[i].push_back(Eigen::Matrix<double, -1, 1>(size));
      grads_ad[i].push_back(Eigen::Matrix<double, -1, 1>(size));
      hess_fd[i].push_back(Eigen::Matrix<double, -1, -1>(size, size));
      hess_ad[i].push_back(Eigen::Matrix<double, -1, -1>(size, size));
    }

  int numels = size + size * (size - 1) / 2;
  Eigen::Matrix<double, -1, 1> x(numels);
  for (int i = 0; i < numels; ++i)
    x(i) = i / 20.0;

  for (size_t i = 0; i < static_cast<size_t>(size); ++i)
    for (size_t j = 0; j < static_cast<size_t>(size); ++j) {
      stan::math::hessian(functowns[i][j], x, evals_ad(i,j), grads_ad[i][j], hess_ad[i][j]);
      stan::math::finite_diff_hessian(functowns[i][j], x,
                                       evals_fd(i,j), grads_fd[i][j], 
                                       hess_fd[i][j]);
      for (int m = 0; m < numels; ++m) 
        for (int n = 0; n < numels; ++n)
          EXPECT_NEAR(hess_fd[i][j](m,n), hess_ad[i][j](m,n), 1e-08);
      EXPECT_FLOAT_EQ(evals_fd(i, j), evals_ad(i, j));
    }
}

void test_grad_hessians(int size) {
  std::vector<std::vector<chol_functor> > functowns;
  std::vector<std::vector<Eigen::Matrix<double, -1, -1> > > hess_ad;
  std::vector<std::vector<Eigen::Matrix<double, -1, -1> > > hess_fd;
  std::vector<std::vector<std::vector<Eigen::Matrix<double, -1, -1> > > > grad_hess_ad;
  std::vector<std::vector<std::vector<Eigen::Matrix<double, -1, -1> > > > grad_hess_fd;
  Eigen::Matrix<double, -1, -1> evals_ad(size,size);
  Eigen::Matrix<double, -1, -1> evals_fd(size,size);
  functowns.resize(size);
  grad_hess_ad.resize(size);
  grad_hess_fd.resize(size);
  hess_ad.resize(size);
  hess_fd.resize(size);

  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) 
      functowns[i].push_back(chol_functor(i, j, size));
    grad_hess_fd[i].resize(size);
    grad_hess_ad[i].resize(size);
    hess_fd[i].resize(size);
    hess_ad[i].resize(size);
  }

  int numels = size + size * (size - 1) / 2;
  Eigen::Matrix<double, -1, 1> x(numels);
  for (int i = 0; i < numels; ++i)
    x(i) = i / 10.0;

  for (size_t i = 0; i < static_cast<size_t>(size); ++i)
    for (size_t j = 0; j < static_cast<size_t>(size); ++j) {
      stan::math::grad_hessian(functowns[i][j], x, evals_ad(i,j), hess_ad[i][j], grad_hess_ad[i][j]);
      stan::math::finite_diff_grad_hessian(functowns[i][j], x,
                                           evals_fd(i,j), hess_fd[i][j], 
                                           grad_hess_fd[i][j]);
      for (size_t k = 0; k < static_cast<size_t>(numels); ++k)
        for (int m = 0; m < numels; ++m) 
          for (int n = 0; n < numels; ++n)
            EXPECT_NEAR(grad_hess_fd[i][j][k](m,n), grad_hess_ad[i][j][k](m,n), 1e-08);
      EXPECT_FLOAT_EQ(evals_fd(i, j), evals_ad(i, j));
    }
}

TEST(AgradMixMatrixCholeskyDecompose, exception_mat_fv) {
  stan::math::matrix_fv m;
  
  m.resize(2,2);
  m << 1.0, 2.0, 
    2.0, 3.0;
  EXPECT_THROW(stan::math::cholesky_decompose(m),std::domain_error);

  m.resize(0, 0);
  EXPECT_NO_THROW(stan::math::cholesky_decompose(m));
  
  m.resize(2, 3);
  EXPECT_THROW(stan::math::cholesky_decompose(m), std::invalid_argument);

  // not symmetric
  m.resize(2,2);
  m << 1.0, 2.0,
    3.0, 4.0;
  EXPECT_THROW(stan::math::cholesky_decompose(m), std::domain_error);
}

TEST(AgradMixMatrixCholeskyDecompose, exception_mat_ffv) {
  stan::math::matrix_ffv m;
  
  m.resize(2,2);
  m << 1.0, 2.0, 
    2.0, 3.0;
  EXPECT_THROW(stan::math::cholesky_decompose(m),std::domain_error);

  m.resize(0, 0);
  EXPECT_NO_THROW(stan::math::cholesky_decompose(m));
  
  m.resize(2, 3);
  EXPECT_THROW(stan::math::cholesky_decompose(m), std::invalid_argument);

  // not symmetric
  m.resize(2,2);
  m << 1.0, 2.0,
    3.0, 4.0;
  EXPECT_THROW(stan::math::cholesky_decompose(m), std::domain_error);
}

TEST(AgradMixMatrixCholeskyDecompose, mat_2nd_deriv) {
  test_hessians(3);
}

TEST(AgradMixMatrixCholeskyDecompose, mat_3rd_deriv) {
  test_grad_hessians(3);
}
