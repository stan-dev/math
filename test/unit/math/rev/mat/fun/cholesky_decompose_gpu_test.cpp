#include <stan/math/rev/mat.hpp>
#include <boost/random/mersenne_twister.hpp>

template<typename T_x>
std::vector<T_x> fill_vec(Eigen::Matrix<T_x, -1, 1> inp) {
  std::vector<T_x> ret_vec;
  ret_vec.reserve(inp.rows());
  for (int i = 0; i < inp.rows(); ++i)
    ret_vec.push_back(inp(i));
  return ret_vec;
}

template<typename T>
Eigen::Matrix<T, -1, -1> create_mat(Eigen::VectorXd inp, 
                                    T alpha,
                                    T len, 
                                    T jitter) {
  std::vector<double> test_inp = fill_vec(inp);
  Eigen::Matrix<T, -1, -1> test_mat_dense = stan::math::cov_exp_quad(test_inp, alpha, len);
  for (int i = 0; i < inp.rows(); ++i)
    test_mat_dense(i,i) = test_mat_dense(i,i) + jitter;
  return test_mat_dense;
}

struct gp_chol {
  Eigen::VectorXd inp, mean, y;
  gp_chol(Eigen::VectorXd inp_, Eigen::VectorXd mean_, Eigen::VectorXd y_) : inp(inp_), mean(mean_), y(y_) { }
  template <typename T>
  T operator()(Eigen::Matrix<T, -1, 1> x) const {
    Eigen::Matrix<T, -1, -1> x_c = create_mat(inp, x[0], x[1], x[2]);
    Eigen::Matrix<T, -1, -1> L = stan::math::cholesky_decompose_gpu(x_c);
    T lp = stan::math::multi_normal_cholesky_lpdf(y, mean, L);
    return lp;
  }
};

struct chol_functor {
  int i, j, K;
  chol_functor(int i_, int j_, int K_) : i(i_), j(j_), K(K_) { }
  template <typename T>
  T operator()(Eigen::Matrix<T, -1, 1> x) const {
    using stan::math::cov_matrix_constrain;
    using stan::math::cholesky_decompose_gpu;
    T lp(0.0);
    Eigen::Matrix<T, -1, -1> x_c = cov_matrix_constrain(x, K, lp);
    Eigen::Matrix<T, -1, -1> L = cholesky_decompose_gpu(x_c);
    lp += L(i, j);
    return lp;
  }
};

struct chol_functor_mult_scal {
  int K;
  Eigen::VectorXd vec;
  chol_functor_mult_scal(int K_, Eigen::VectorXd vec_) : K(K_), vec(vec_) { }
  template <typename T>
  T operator()(Eigen::Matrix<T, -1, 1> x) const {
    using stan::math::cov_matrix_constrain;
    using stan::math::cholesky_decompose_gpu;
    using stan::math::multiply;
    using stan::math::transpose;
    T lp(0.0);
    Eigen::Matrix<T, -1, -1> x_c = cov_matrix_constrain(x, K, lp);
    Eigen::Matrix<T, -1, -1> L = cholesky_decompose_gpu(x_c);
    lp += multiply(transpose(vec), multiply(L, vec));
    return lp;
  }
};


struct chol_functor_2 {
  int K;
  chol_functor_2(int K_) : K(K_) { }
  template <typename T>
  T operator()(Eigen::Matrix<T, -1, 1> x) const {
    using stan::math::cov_matrix_constrain;
    using stan::math::cholesky_decompose_gpu;
    using stan::math::multi_normal_cholesky_log;
    T lp(0.0);
    Eigen::Matrix<T, -1, -1> x_c = cov_matrix_constrain(x, K, lp);
    Eigen::Matrix<T, -1, -1> L = cholesky_decompose_gpu(x_c);
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
    using stan::math::cholesky_decompose_gpu;
    Eigen::Matrix<T, -1, -1> x_c(K, K);
    int pos = 0;
    for (int n = 0; n < K; ++n) 
      for (int m = 0; m < K; ++m) {
        x_c(m,n) = x(pos++);
        x_c(n,m) = x_c(m,n);
      }
    Eigen::Matrix<T, -1, -1> L = cholesky_decompose_gpu(x_c);
    return L(i, j);
  }
};


void test_gradients(int size, double prec) {
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
      x(i) = i % 10 / 100.0 ;

  for (size_t i = 0; i < static_cast<size_t>(size); ++i) {
    for (size_t j = 0; j < static_cast<size_t>(size); ++j) {
      stan::math::gradient(functors[i][j], x, evals_ad(i,j), grads_ad[i][j]);
      stan::math::finite_diff_gradient(functors[i][j], x,
                                       evals_fd(i,j), grads_fd[i][j]);
    
    }
  }
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
  stan::math::finite_diff_gradient(functown, x,evals_fd, grads_fd);

  return grads_ad.sum();
}


int main() {
  test_gradient(1500, 1e-10);
  test_gradients(1009, 1e-10);
  return 1;
}

