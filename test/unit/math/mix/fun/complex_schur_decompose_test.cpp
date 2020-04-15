#include <test/unit/math/test_ad.hpp>
#include <stdexcept>

TEST(mathMixFun, complexSchurDecomposeT) {
  auto f = [](const auto& x) {
    using stan::math::complex_schur_decompose_t;
    return complex_schur_decompose_t(x);
  };
  for (const auto& x : stan::test::square_test_matrices(0, 2)) {
    stan::test::expect_ad(f, x);
  }

  Eigen::MatrixXd a32(3, 2);
  a32 << 3, -5, 7, -7.2, 9.1, -6.3;
  EXPECT_THROW(f(a32), std::invalid_argument);
}

TEST(mathMixFun, complexSchurDecomposeU) {
  auto f = [](const auto& x) {
    using stan::math::complex_schur_decompose_u;
    return complex_schur_decompose_u(x);
  };
  for (const auto& x : stan::test::square_test_matrices(0, 2)) {
    stan::test::expect_ad(f, x);
  }

  Eigen::MatrixXd a32(3, 2);
  a32 << 3, -5, 7, -7.2, 9.1, -6.3;
  EXPECT_THROW(f(a32), std::invalid_argument);
}

template <typename V>
void test_complex_schur_decompose(const Eigen::MatrixXd& x) {
  using stan::math::complex_schur_decompose_t;
  using stan::math::complex_schur_decompose_u;
  using stan::math::value_of_rec;
  Eigen::Matrix<V, -1, -1> X(x.rows(), x.cols());
  for (int i = 0; i < x.size(); ++i)
    X(i) = x(i);
  auto T = complex_schur_decompose_t(X);
  auto U = complex_schur_decompose_u(X);
  auto X2 = U * T * U.adjoint();
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < x.rows(); ++i) {
      EXPECT_FLOAT_EQ(x(i, j), value_of_rec(X2(i, j).real()));
    }
  }
}

// TODO(carpenter): complex_schur_decompose only works for real types,
// returning complex types;  the functions work for double types, but
// not autodiff complex types;  what's missing is compound arithmetic
// for matrix base type lhs and complex rhs

template <typename V>
void test_complex_schur_decompose_complex(const Eigen::MatrixXd& x) {
  using stan::math::complex_schur_decompose_t;
  using stan::math::complex_schur_decompose_u;
  using stan::math::value_of_rec;
  Eigen::Matrix<std::complex<V>, -1, -1> X(x.rows(), x.cols());
  for (int i = 0; i < x.size(); ++i)
    X(i) = std::complex<double>(x(i), i);
  auto T = complex_schur_decompose_t(X);
  auto U = complex_schur_decompose_u(X);
  auto X2 = U * T * U.adjoint();
  for (int j = 0; j < x.cols(); ++j) {
    for (int i = 0; i < x.rows(); ++i) {
      EXPECT_NEAR(value_of_rec(X(i, j).real()), value_of_rec(X2(i, j).real()),
                  1e-8);
      EXPECT_NEAR(value_of_rec(X(i, j)).imag(), value_of_rec(X2(i, j).imag()),
                  1e-8);
    }
  }
}

TEST(mathMixFun, complexSchurDecompose) {
  using d_t = double;
  using v_t = stan::math::var;
  using fd_t = stan::math::fvar<d_t>;
  using ffd_t = stan::math::fvar<fd_t>;
  using fv_t = stan::math::fvar<v_t>;
  using ffv_t = stan::math::fvar<fv_t>;
  for (const auto& x : stan::test::square_test_matrices(1, 3)) {
    test_complex_schur_decompose<d_t>(x);
    test_complex_schur_decompose<v_t>(x);
    test_complex_schur_decompose<fd_t>(x);
    test_complex_schur_decompose<ffd_t>(x);
    test_complex_schur_decompose<fv_t>(x);
    test_complex_schur_decompose<ffv_t>(x);
  }
  for (const auto& x : stan::test::square_test_matrices(1, 3)) {
    test_complex_schur_decompose_complex<d_t>(x);
    // TODO(carpenter): // require expanding compound arithmetic
    //                  // to more general Eigen types for *=, -=, etc.
    // test_complex_schur_decompose_complex<v_t>(x);
    // test_complex_schur_decompose_complex<fd_t>(x);
    // test_complex_schur_decompose_complex<ffd_t>(x);
    // test_complex_schur_decompose_complex<fv_t>(x);
    // test_complex_schur_decompose_complex<ffv_t>(x);
  }
}
