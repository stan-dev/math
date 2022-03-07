#include <test/unit/math/test_ad.hpp>

template <typename T>
Eigen::Matrix<T, -1, 1>
flatten(const Eigen::Matrix<std::complex<T>, -1, 1>& x) {
  Eigen::Matrix<T, -1, 1> xf(x.size() / 2);
  for (int i = 0; i < x.size(); ++i) {
    xf[i / 2] = std::real(xf[i]);
    xf[i / 2 + 1] = std::imag(xf[i]);
  }
  return xf;
}

template <typename T>
Eigen::Matrix<T, -1, 1> to_real_vec(const Eigen::Matrix<std::complex<T>, -1, 1>& x) {
  Eigen::Matrix<T, -1, 1> y(x.size() * 2);
  for (int i = 0; i < y.size(); i += 2) {
    y[i] = real(x[i / 2]);
    y[i + 1] = imag(x[i / 2]);
  }
  return y;
}

template <typename T>
Eigen::Matrix<std::complex<T>, -1, 1> to_complex_vec(const Eigen::Matrix<T, -1, 1>& x) {
  Eigen::Matrix<std::complex<T>, -1, 1> y(x.size() / 2);
  for (int i = 0; i < y.size(); ++i)
    y[i] = { x[2 * i], x[2 * i + 1] };
  return y;
}

struct foo {
  template <typename T>
  Eigen::Matrix<T, -1, 1> operator()(const Eigen::Matrix<T, -1, 1>& x) const {
    auto xc = to_complex_vec(x);
    auto y = stan::math::fft(xc);
    return to_real_vec(y);
  }
};    




TEST(mathMixFun, fft) {
  typedef Eigen::Matrix<std::complex<double>, -1, 1> cvec_t;
  typedef std::complex<double> c_t;
  auto f = [](const auto& x) {
    using stan::math::fft;
    return fft(x);
  };
  
  cvec_t x0(0);
  stan::test::expect_ad(f, x0);

  cvec_t x1(1);
  x1[0] = { 1, 2 };
  stan::test::expect_ad(f, x1);

  cvec_t x2(2);
  x2[0] = { 1, 2};
  x2[1] = {-1.3, 2.9};
  
  // stan::test::expect_ad(f, x2);

  Eigen::VectorXd x(6);
  x << 1, -1.3, 2.9, 14.7, -12.9, -4.8;

  Eigen::VectorXd gx;
  Eigen::MatrixXd J;
  foo g;
  stan::math::jacobian(g, x, gx, J);
  std::cout << "J = " << std::endl << J << std::endl;
  
}

