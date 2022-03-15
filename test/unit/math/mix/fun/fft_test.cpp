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
  int n = 0; int i = 0;
  auto f = [&](const auto& x) {
    using stan::math::fft;
    return i ? real(fft(x)(n)) : imag(fft(x)(n));
  };

  cvec_t x0(0);
  EXPECT_EQ(0, stan::math::fft(x0).size());  // no AD to test

  cvec_t x1(1);
  x1[0] = { 1, 2 };
  for (i = 0; i <= 1; ++i)
    for (n = 0; n <= 0; ++n)
      stan::test::expect_ad(f, x1);

  cvec_t x2(2);
  x2[0] = {1, 2};
  x2[1] = {-1.3, 2.9};
  for (i = 0; i <= 1; ++i)
    for (n = 0; n <= 1; ++n)
      stan::test::expect_ad(f, x2);

  Eigen::VectorXcd x3(3);
  x3[0] = {1, -1.3};
  x3[1] = {2.9, 14.7};
  x3[2] = {-12.9, -4.8};
  for (i = 0; i <= 1; ++i)
    for (n = 0; n <= 2; ++n)
      stan::test::expect_ad(f, x3);

  Eigen::VectorXcd x4(4);
  x4[0] = {1, -1.3};
  x4[1] = {-2.9, 14.7};
  x4[2] = {-12.9, -4.8};
  x4[3] = {8.398, 9.387};
  for (i = 0; i <= 1; ++i)
    for (n = 0; n <= 3; ++n)
      stan::test::expect_ad(f, x4);
  
}

