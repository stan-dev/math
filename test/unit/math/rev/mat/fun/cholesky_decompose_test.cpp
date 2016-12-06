#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <boost/random/mersenne_twister.hpp>

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

struct chol_functor_2 {
  int K;
  chol_functor_2(int K_) : K(K_) { }
  template <typename T>
  T operator()(Eigen::Matrix<T, -1, 1> x) const {
    using stan::math::cov_matrix_constrain;
    using stan::math::cholesky_decompose;
    using stan::math::multi_normal_cholesky_log;
    T lp(0.0);
    Eigen::Matrix<T, -1, -1> x_c = cov_matrix_constrain(x, K, lp);
    Eigen::Matrix<T, -1, -1> L = cholesky_decompose(x_c);
    Eigen::Matrix<double, -1, 1> vec(K);
    Eigen::Matrix<double, -1, 1> mu(K);
    vec.setZero();
    mu.setOnes();
    lp += multi_normal_cholesky_log(vec, mu, L);
    return lp;
  }
};

struct chol_functor_simple {
  int i, j, K;
  chol_functor_simple(int i_, int j_, int K_) : i(i_), j(j_), K(K_) { }
  template <typename T>
  T operator()(Eigen::Matrix<T, -1, 1> x) const {
    using stan::math::cholesky_decompose;
    Eigen::Matrix<T, -1, -1> x_c(K, K);
    int pos = 0;
    for (int n = 0; n < K; ++n) 
      for (int m = 0; m < K; ++m) {
        x_c(m,n) = x(pos++);
        x_c(n,m) = x_c(m,n);
      }
    Eigen::Matrix<T, -1, -1> L = cholesky_decompose(x_c);
    return L(i, j);
  }
};


void test_gradients(int size) {
  std::vector<std::vector<chol_functor> > functors;
  std::vector<std::vector<Eigen::Matrix<double, -1, 1> > > grads_ad;
  std::vector<std::vector<Eigen::Matrix<double, -1, 1> > > grads_fd;
  Eigen::Matrix<double, -1, -1> evals_ad(size,size);
  Eigen::Matrix<double, -1, -1> evals_fd(size,size);
  functors.resize(size);
  grads_ad.resize(size);
  grads_fd.resize(size);

  for (int i = 0; i < size; ++i)
    for (int j = 0; j < size; ++j) {
      functors[i].push_back(chol_functor(i, j, size));
      grads_fd[i].push_back(Eigen::Matrix<double, -1, 1>(size));
      grads_ad[i].push_back(Eigen::Matrix<double, -1, 1>(size));
    }

  int numels = size + size * (size - 1) / 2;
  Eigen::Matrix<double, -1, 1> x(numels);
  for (int i = 0; i < numels; ++i)
    x(i) = i / 10.0;

  for (size_t i = 0; i < static_cast<size_t>(size); ++i) {
    for (size_t j = 0; j < static_cast<size_t>(size); ++j) {
      stan::math::gradient(functors[i][j], x, evals_ad(i,j), grads_ad[i][j]);
      stan::math::finite_diff_gradient(functors[i][j], x,
                                       evals_fd(i,j), grads_fd[i][j]);
    
      for (int k = 0; k < numels; ++k) 
        EXPECT_NEAR(grads_fd[i][j](k), grads_ad[i][j](k), 1e-10);
      EXPECT_FLOAT_EQ(evals_fd(i, j), evals_ad(i, j));
    }
  }
}

void test_gradients_simple(int size) {
  std::vector<std::vector<chol_functor_simple> > functors;
  std::vector<std::vector<Eigen::Matrix<double, -1, 1> > > grads_ad;
  std::vector<std::vector<Eigen::Matrix<double, -1, 1> > > grads_fd;
  Eigen::Matrix<double, -1, -1> evals_ad(size,size);
  Eigen::Matrix<double, -1, -1> evals_fd(size,size);
  functors.resize(size);
  grads_ad.resize(size);
  grads_fd.resize(size);

  for (int i = 0; i < size; ++i)
    for (int j = 0; j < size; ++j) {
      functors[i].push_back(chol_functor_simple(i, j, size));
      grads_fd[i].push_back(Eigen::Matrix<double, -1, 1>(size));
      grads_ad[i].push_back(Eigen::Matrix<double, -1, 1>(size));
    }

  stan::math::welford_covar_estimator estimator(size);
  
  boost::random::mt19937 rng;
  for (int i = 0; i < 1000; ++i) {
    Eigen::VectorXd q(size);
    for (int j = 0; j < size; ++j)
      q(j) = stan::math::normal_rng(0.0,1.0,rng);
    estimator.add_sample(q);
  }

  Eigen::MatrixXd covar(size, size);
  estimator.sample_covariance(covar);

  Eigen::Matrix<double, -1, 1> x(size * size);
  int pos = 0;
  for (int j = 0; j < size; ++j)
    for (int i = 0; i < size; ++i)
      x(pos++) = covar(i, j);

  for (size_t j = 0; j < static_cast<size_t>(size); ++j) {
    for (size_t i = j; i < static_cast<size_t>(size); ++i) {
      stan::math::gradient(functors[i][j], x, evals_ad(i,j), grads_ad[i][j]);
      stan::math::finite_diff_gradient(functors[i][j], x,
                                       evals_fd(i,j), grads_fd[i][j]);
    
      for (int k = 0; k < size; ++k) 
        EXPECT_NEAR(grads_fd[i][j](k), grads_ad[i][j](k), 1e-10);
      EXPECT_FLOAT_EQ(evals_fd(i, j), evals_ad(i, j));
    }
  }
}

double test_gradient(int size) {
  chol_functor_2 functown(size);
  Eigen::Matrix<double, -1, 1> grads_ad;
  Eigen::Matrix<double, -1, 1> grads_fd;
  double evals_ad;
  double evals_fd;

  int numels = size + size * (size - 1) / 2;
  Eigen::Matrix<double, -1, 1> x(numels);
  for (int i = 0; i < numels; ++i)
    x(i) = i / 100.0;

  stan::math::gradient(functown, x, evals_ad, grads_ad);
  stan::math::finite_diff_gradient(functown, x,evals_fd, grads_fd);

  for (int k = 0; k < numels; ++k) 
    EXPECT_NEAR(grads_fd(k), grads_ad(k), 1e-09);
  EXPECT_FLOAT_EQ(evals_fd, evals_ad);
  return grads_ad.sum();
}

TEST(AgradRevMatrix,mat_cholesky) {
  using stan::math::matrix_v;
  using stan::math::transpose;
  using stan::math::cholesky_decompose;
  using stan::math::singular_values;

  // symmetric
  matrix_v X(2,2);
  AVAR a = 3.0;
  AVAR b = -1.0;
  AVAR c = -1.0;
  AVAR d = 1.0;
  X << a, b, 
    c, d;

  matrix_v L = cholesky_decompose(X);

  matrix_v LL_trans = multiply(L,transpose(L));
  EXPECT_FLOAT_EQ(a.val(),LL_trans(0,0).val());
  EXPECT_FLOAT_EQ(b.val(),LL_trans(0,1).val());
  EXPECT_FLOAT_EQ(c.val(),LL_trans(1,0).val());
  EXPECT_FLOAT_EQ(d.val(),LL_trans(1,1).val());

  EXPECT_NO_THROW(singular_values(X));
}

TEST(AgradRevMatrix, exception_mat_cholesky) {
  stan::math::matrix_v m;
  
  // not positive definite
  m.resize(2,2);
  m << 1.0, 2.0, 
    2.0, 3.0;
  EXPECT_THROW(stan::math::cholesky_decompose(m),std::domain_error);

  // zero size
  m.resize(0, 0);
  EXPECT_NO_THROW(stan::math::cholesky_decompose(m));
  
  // not square
  m.resize(2, 3);
  EXPECT_THROW(stan::math::cholesky_decompose(m), std::invalid_argument);

  // not symmetric
  m.resize(2,2);
  m << 1.0, 2.0,
    3.0, 4.0;
  EXPECT_THROW(stan::math::cholesky_decompose(m), std::domain_error);
}

TEST(AgradRevMatrix, mat_cholesky_1st_deriv) {
  test_gradients(9);
  test_gradients_simple(10);
  test_gradient(50);
}

TEST(AgradRevMatrix, check_varis_on_stack) {
  stan::math::matrix_v X(2,2);
  X << 3, -1, -1, 1;

  test::check_varis_on_stack(stan::math::cholesky_decompose(X));
}
