#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/util.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <vector>

template <typename T_x>
std::vector<T_x> fill_vec(Eigen::Matrix<T_x, -1, 1> inp) {
  std::vector<T_x> ret_vec;
  ret_vec.reserve(inp.rows());
  for (int i = 0; i < inp.rows(); ++i)
    ret_vec.push_back(inp(i));
  return ret_vec;
}

template <typename T>
Eigen::Matrix<T, -1, -1> create_mat(Eigen::VectorXd inp, T alpha, T len,
                                    T jitter) {
  std::vector<double> test_inp = fill_vec(inp);
  Eigen::Matrix<T, -1, -1> test_mat_dense
      = stan::math::cov_exp_quad(test_inp, alpha, len);
  for (int i = 0; i < inp.rows(); ++i)
    test_mat_dense(i, i) = test_mat_dense(i, i) + jitter;
  return test_mat_dense;
}

struct gp_chol {
  Eigen::VectorXd inp, mean, y;
  gp_chol(Eigen::VectorXd inp_, Eigen::VectorXd mean_, Eigen::VectorXd y_)
      : inp(inp_), mean(mean_), y(y_) {}
  template <typename T>
  T operator()(Eigen::Matrix<T, -1, 1> x) const {
    Eigen::Matrix<T, -1, -1> x_c = create_mat(inp, x[0], x[1], x[2]);
    Eigen::Matrix<T, -1, -1> L = stan::math::cholesky_decompose(x_c);
    T lp = stan::math::multi_normal_cholesky_lpdf(y, mean, L);
    return lp;
  }
};

struct chol_functor {
  int i, j, K;
  chol_functor(int i_, int j_, int K_) : i(i_), j(j_), K(K_) {}
  template <typename T>
  T operator()(Eigen::Matrix<T, -1, 1> x) const {
    using stan::math::cholesky_decompose;
    using stan::math::cov_matrix_constrain;
    T lp(0.0);
    Eigen::Matrix<T, -1, -1> x_c = cov_matrix_constrain(x, K, lp);
    Eigen::Matrix<T, -1, -1> L = cholesky_decompose(x_c);
    lp += L(i, j);
    return lp;
  }
};

struct chol_functor_mult_scal {
  int K;
  Eigen::VectorXd vec;
  chol_functor_mult_scal(int K_, Eigen::VectorXd vec_) : K(K_), vec(vec_) {}
  template <typename T>
  T operator()(Eigen::Matrix<T, -1, 1> x) const {
    using stan::math::cholesky_decompose;
    using stan::math::cov_matrix_constrain;
    using stan::math::multiply;
    using stan::math::transpose;
    T lp(0.0);
    Eigen::Matrix<T, -1, -1> x_c = cov_matrix_constrain(x, K, lp);
    Eigen::Matrix<T, -1, -1> L = cholesky_decompose(x_c);
    lp += multiply(transpose(vec), multiply(L, vec));
    return lp;
  }
};

struct chol_functor_2 {
  int K;
  explicit chol_functor_2(int K_) : K(K_) {}
  template <typename T>
  T operator()(Eigen::Matrix<T, -1, 1> x) const {
    using stan::math::cholesky_decompose;
    using stan::math::cov_matrix_constrain;
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
  chol_functor_simple(int i_, int j_, int K_) : i(i_), j(j_), K(K_) {}
  template <typename T>
  T operator()(Eigen::Matrix<T, -1, 1> x) const {
    using stan::math::cholesky_decompose;
    Eigen::Matrix<T, -1, -1> x_c(K, K);
    int pos = 0;
    for (int n = 0; n < K; ++n)
      for (int m = 0; m < K; ++m) {
        x_c(m, n) = x(pos++);
        x_c(n, m) = x_c(m, n);
      }
    Eigen::Matrix<T, -1, -1> L = cholesky_decompose(x_c);
    return L(i, j);
  }
};

struct chol_functor_simple_vec {
  int K;
  Eigen::VectorXd vec;
  chol_functor_simple_vec(int K_, Eigen::VectorXd vec_) : K(K_), vec(vec_) {}
  template <typename T>
  T operator()(Eigen::Matrix<T, -1, 1> x) const {
    using stan::math::cholesky_decompose;
    using stan::math::multiply;
    using stan::math::transpose;
    Eigen::Matrix<T, -1, -1> x_c(K, K);
    int pos = 0;
    for (int n = 0; n < K; ++n)
      for (int m = 0; m < K; ++m) {
        x_c(m, n) = x(pos++);
        x_c(n, m) = x_c(m, n);
      }
    Eigen::Matrix<T, -1, -1> L = cholesky_decompose(x_c);
    T lp = multiply(transpose(vec), multiply(L, vec));
    return lp;
  }
};

void test_gradients(int size, double prec) {
  std::vector<std::vector<chol_functor> > functors;
  std::vector<std::vector<Eigen::Matrix<double, -1, 1> > > grads_ad;
  std::vector<std::vector<Eigen::Matrix<double, -1, 1> > > grads_fd;
  Eigen::Matrix<double, -1, -1> evals_ad(size, size);
  Eigen::Matrix<double, -1, -1> evals_fd(size, size);
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
    x(i) = i % 10 / 100.0;

  for (size_t i = 0; i < static_cast<size_t>(size); ++i) {
    for (size_t j = 0; j < static_cast<size_t>(size); ++j) {
      stan::math::gradient(functors[i][j], x, evals_ad(i, j), grads_ad[i][j]);
      stan::math::finite_diff_gradient(functors[i][j], x, evals_fd(i, j),
                                       grads_fd[i][j]);

      for (int k = 0; k < numels; ++k)
        EXPECT_NEAR(grads_fd[i][j](k), grads_ad[i][j](k), prec);
      EXPECT_FLOAT_EQ(evals_fd(i, j), evals_ad(i, j));
    }
  }
}

void test_gradients_simple(int size, double prec) {
  std::vector<std::vector<chol_functor_simple> > functors;
  std::vector<std::vector<Eigen::Matrix<double, -1, 1> > > grads_ad;
  std::vector<std::vector<Eigen::Matrix<double, -1, 1> > > grads_fd;
  Eigen::Matrix<double, -1, -1> evals_ad(size, size);
  Eigen::Matrix<double, -1, -1> evals_fd(size, size);
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
      q(j) = stan::math::normal_rng(0.0, 1.0, rng);
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
      stan::math::gradient(functors[i][j], x, evals_ad(i, j), grads_ad[i][j]);
      stan::math::finite_diff_gradient(functors[i][j], x, evals_fd(i, j),
                                       grads_fd[i][j]);

      for (int k = 0; k < size; ++k)
        EXPECT_NEAR(grads_fd[i][j](k), grads_ad[i][j](k), prec);
      EXPECT_FLOAT_EQ(evals_fd(i, j), evals_ad(i, j));
    }
  }
}

void test_gp_grad(int mat_size, double prec) {
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::var;

  Eigen::VectorXd test_vec(mat_size);
  Eigen::VectorXd mean_vec(mat_size);
  Eigen::VectorXd draw_vec(mat_size);
  Eigen::VectorXd test_vals(3);
  test_vals[0] = 1;
  test_vals[1] = 1.5;
  test_vals[2] = 10;

  boost::random::mt19937 rng(2);

  for (int i = 0; i < mat_size; ++i) {
    test_vec(i) = stan::math::normal_rng(0.0, 0.1, rng);
    mean_vec(i) = 0;
    draw_vec(i) = stan::math::normal_rng(0.0, 0.1, rng);
  }

  Eigen::MatrixXd cov_mat
      = create_mat(test_vec, test_vals[0], test_vals[1], test_vals[2]);
  Eigen::MatrixXd chol_cov = cov_mat.llt().matrixL();
  Eigen::VectorXd y_vec = chol_cov * draw_vec;

  gp_chol gp_fun(test_vec, mean_vec, y_vec);
  double val_ad;
  Eigen::VectorXd grad_ad;
  stan::math::gradient(gp_fun, test_vals, val_ad, grad_ad);

  VectorXd grad_fd;
  double val_fd;

  stan::math::finite_diff_gradient(gp_fun, test_vals, val_fd, grad_fd);
  EXPECT_NEAR(val_fd, val_ad, 1e-10);
  for (int i = 0; i < grad_ad.size(); ++i) {
    EXPECT_NEAR(grad_fd(i), grad_ad(i), prec);
  }
}

void test_chol_mult(int mat_size, double prec) {
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::math::var;

  int vec_size = mat_size * (mat_size + 1) / 2;
  Eigen::VectorXd test_vec(mat_size);
  Eigen::VectorXd test_vals(vec_size);

  boost::random::mt19937 rng(2);

  for (int i = 0; i < test_vals.size(); ++i) {
    if (i < test_vec.size()) {
      test_vec(i) = stan::math::normal_rng(0.0, 0.1, rng);
    }
    test_vals(i) = i % 10 / 100.0;
  }

  chol_functor_mult_scal mult_fun(mat_size, test_vec);
  double val_ad;
  Eigen::VectorXd grad_ad;
  stan::math::gradient(mult_fun, test_vals, val_ad, grad_ad);

  VectorXd grad_fd;
  double val_fd;

  stan::math::finite_diff_gradient(mult_fun, test_vals, val_fd, grad_fd);
  EXPECT_NEAR(val_fd, val_ad, 1e-10);
  for (int i = 0; i < grad_ad.size(); ++i) {
    EXPECT_NEAR(grad_fd(i), grad_ad(i), prec);
  }
}

void test_simple_vec_mult(int size, double prec) {
  Eigen::VectorXd test_vec(size);
  boost::random::mt19937 rng(2);

  for (int i = 0; i < test_vec.size(); ++i)
    test_vec(i) = stan::math::normal_rng(0.0, 0.1, rng);

  chol_functor_simple_vec f(size, test_vec);

  stan::math::welford_covar_estimator estimator(size);

  for (int i = 0; i < 1000; ++i) {
    Eigen::VectorXd q(size);
    for (int j = 0; j < size; ++j)
      q(j) = stan::math::normal_rng(0.0, 1.0, rng);
    estimator.add_sample(q);
  }

  Eigen::MatrixXd covar(size, size);
  estimator.sample_covariance(covar);

  Eigen::Matrix<double, -1, 1> x(size * size);
  int pos = 0;
  for (int j = 0; j < size; ++j)
    for (int i = 0; i < size; ++i)
      x(pos++) = covar(i, j);

  double eval_ad;
  Eigen::VectorXd grad_ad;
  stan::math::gradient(f, x, eval_ad, grad_ad);
  double eval_fd;
  Eigen::VectorXd grad_fd;
  stan::math::finite_diff_gradient(f, x, eval_fd, grad_fd);

  EXPECT_FLOAT_EQ(eval_fd, eval_ad);
  for (int k = 0; k < grad_fd.size(); ++k)
    EXPECT_NEAR(grad_fd(k), grad_ad(k), prec) << " for k=" << k;
}

double test_gradient(int size, double prec) {
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
  stan::math::finite_diff_gradient(functown, x, evals_fd, grads_fd);

  for (int k = 0; k < numels; ++k)
    EXPECT_NEAR(grads_fd(k), grads_ad(k), prec) << " for k=" << k;
  EXPECT_FLOAT_EQ(evals_fd, evals_ad);
  return grads_ad.sum();
}

TEST(AgradRevMatrix, mat_cholesky) {
  using stan::math::cholesky_decompose;
  using stan::math::matrix_v;
  using stan::math::singular_values;
  using stan::math::transpose;

  // symmetric
  matrix_v X(2, 2);
  stan::math::var a = 3.0;
  stan::math::var b = -1.0;
  stan::math::var c = -1.0;
  stan::math::var d = 1.0;
  X << a, b, c, d;

  matrix_v L = cholesky_decompose(X);

  matrix_v LL_trans = multiply(L, transpose(L));
  EXPECT_FLOAT_EQ(a.val(), LL_trans(0, 0).val());
  EXPECT_FLOAT_EQ(b.val(), LL_trans(0, 1).val());
  EXPECT_FLOAT_EQ(c.val(), LL_trans(1, 0).val());
  EXPECT_FLOAT_EQ(d.val(), LL_trans(1, 1).val());

  EXPECT_NO_THROW(singular_values(X));
}

TEST(AgradRevMatrix, exception_mat_cholesky) {
  stan::math::matrix_v m;

  // not positive definite
  m.resize(2, 2);
  m << 1.0, 2.0, 2.0, 3.0;
  EXPECT_THROW(stan::math::cholesky_decompose(m), std::domain_error);

  // zero size
  m.resize(0, 0);
  EXPECT_NO_THROW(stan::math::cholesky_decompose(m));

  // not square
  m.resize(2, 3);
  EXPECT_THROW(stan::math::cholesky_decompose(m), std::invalid_argument);

  // not symmetric
  m.resize(2, 2);
  m << 1.0, 2.0, 3.0, 4.0;
  EXPECT_THROW(stan::math::cholesky_decompose(m), std::domain_error);
}

TEST(AgradRevMatrix, exception_varmat_cholesky) {
  stan::math::matrix_d m;

  // not positive definite
  m.resize(2, 2);
  m << 1.0, 2.0, 2.0, 3.0;
  stan::math::var_value<stan::math::matrix_d> mv1(m);
  EXPECT_THROW(stan::math::cholesky_decompose(mv1), std::domain_error);

  // zero size
  m.resize(0, 0);
  stan::math::var_value<stan::math::matrix_d> mv2(m);
  EXPECT_NO_THROW(stan::math::cholesky_decompose(mv2));

  // not square
  m.resize(2, 3);
  stan::math::var_value<stan::math::matrix_d> mv3(m);
  EXPECT_THROW(stan::math::cholesky_decompose(mv3), std::invalid_argument);

  // not symmetric
  m.resize(2, 2);
  m << 1.0, 2.0, 3.0, 4.0;
  stan::math::var_value<stan::math::matrix_d> mv4(m);
  EXPECT_THROW(stan::math::cholesky_decompose(mv4), std::domain_error);
}

TEST(AgradRevMatrix, mat_cholesky_1st_deriv_small) {
  test_gradients(9, 1e-10);
  test_gradients_simple(10, 1e-10);
  test_gradient(15, 1e-10);
  test_gp_grad(20, 1e-10);
}

TEST(AgradRevMatrix, check_varis_on_stack_small) {
  stan::math::matrix_v X(2, 2);
  X << 3, -1, -1, 1;

  test::check_varis_on_stack(stan::math::cholesky_decompose(X));
}

TEST(AgradRevMatrix, mat_cholesky_1st_deriv_large_gradients) {
  test_gradient(36, 1e-08);
  test_gp_grad(100, 1e-08);
  test_gp_grad(1000, 1e-08);
  test_chol_mult(37, 1e-08);
  test_simple_vec_mult(45, 1e-08);
}

TEST(AgradRevMatrix, cholesky_replicated_input) {
  using stan::math::var;

  auto f = [](int size, const auto& y) {
    auto m = stan::math::diag_matrix(stan::math::rep_vector(y, size));
    auto L = stan::math::cholesky_decompose(m);
    return stan::math::sum(L);
  };

  double ydbl = 1.5;
  double dx = 1e-5;
  var y = ydbl;
  int size = 4;
  var s = f(size, y);
  s.grad();

  double fd_ref = (f(size, ydbl + dx) - f(size, ydbl - dx)) / (2.0 * dx);
  EXPECT_FLOAT_EQ(y.adj(), fd_ref);

  stan::math::set_zero_all_adjoints();
  size = 40;
  s = f(size, y);
  s.grad();

  fd_ref = (f(size, ydbl + dx) - f(size, ydbl - dx)) / (2.0 * dx);
  EXPECT_FLOAT_EQ(y.adj(), fd_ref);
}
