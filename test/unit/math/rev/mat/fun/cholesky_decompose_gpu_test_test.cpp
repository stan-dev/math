#include <stan/math/rev/mat.hpp>
#include <boost/random/mersenne_twister.hpp>


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
 //   std::cout << " INPUT: \n" << x_c << "\n";
    Eigen::Matrix<T, -1, -1> L = cholesky_decompose_gpu(x_c);
 //   std::cout << " OUTPUT: \n" << L << "\n";
    Eigen::Matrix<double, -1, 1> vec(K);
    Eigen::Matrix<double, -1, 1> mu(K);
    vec.setZero();
    mu.setOnes();
    lp += multi_normal_cholesky_log(vec, mu, L);
    return lp;
  }
};


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
  test_gradient(36, 1e-08);
  //test_gp_grad(100, 1e-08);
  //test_gp_grad(1000, 1e-08);
  //test_chol_mult(37, 1e-08);
  //test_simple_vec_mult(45, 1e-08);
  return 1;
}

