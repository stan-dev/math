#include <test/unit/math/test_ad.hpp>

TEST(mathMixFun, rfftroundtrip) {
  using stan::math::inv_rfft2;
  using stan::math::rfft2;
  for (const auto& m_d : stan::test::square_test_matrices(0, 10)) {
    auto fourier_space = rfft2(m_d);
    auto roundtripped = inv_rfft2(fourier_space);
    EXPECT_MATRIX_NEAR(roundtripped, m_d, 1e-8);
  }
}

void expect_rfft(const Eigen::VectorXd& x) {
  for (int m = 0; m < x.rows(); ++m) {
    auto g = [m](const auto& x) {
      using stan::math::rfft;
      return rfft(x)(m);
    };
    stan::test::expect_ad(g, x);
  }
}

TEST(mathMixFun, rfft) {
  using vec_t = Eigen::VectorXd;

  vec_t x0(0);
  expect_rfft(x0);

  vec_t x1(1);
  x1 << 1;
  expect_rfft(x1);

  vec_t x2(2);
  x2 << 1, 2.9;
  expect_rfft(x2);

  vec_t x3(3);
  x3 << 1, 2.9, -14.7;
  expect_rfft(x3);

  vec_t x4(4);
  x4 << 1, 2.9, -14.7, 9.734;
  expect_rfft(x4);
}

void expect_inv_rfft(const Eigen::VectorXcd& x) {
  for (int m = 0; m < x.rows(); ++m) {
    auto g = [m](const auto& x) {
      using stan::math::inv_rfft;
      return inv_rfft(x)(m);
    };
    stan::test::expect_ad(g, x);
  }
}

TEST(mathMixFun, inv_rfft) {
  using cvec_t = Eigen::VectorXcd;
  using c_t = std::complex<double>;

  cvec_t x0(0);
  expect_inv_rfft(x0);

  cvec_t x1(1);
  x1[0] = {-3.247, 0};
  expect_inv_rfft(x1);

  cvec_t x3(3);
  x3[0] = {5, 0};
  x3[1] = {-1, 8.6602540};
  x3[2] = {-1, -8.6602540};
  expect_inv_rfft(x3);

  cvec_t x4(4);
  x4 << c_t(-1.066, -0.), c_t(15.7, +6.834), c_t(-26.334, -0.),
      c_t(15.7, -6.834);
  expect_inv_rfft(x4);
}

void expect_rfft2(const Eigen::MatrixXd& x) {
  for (int n = 0; n < x.cols(); ++n) {
    for (int m = 0; m < x.rows(); ++m) {
      auto g = [m, n](const auto& x) {
        using stan::math::rfft2;
        return rfft2(x)(m, n);
      };
      stan::test::expect_ad(g, x);
    }
  }
}

TEST(mathMixFun, rfft2) {
  Eigen::MatrixXd x(0, 0);
  expect_rfft2(x);

  Eigen::MatrixXd x11(1, 1);
  x11 << -3.9;
  expect_rfft2(x11);

  //   Eigen::MatrixXd x12(1, 2);
  //   x12 << -1.9, 0.2;
  //   expect_rfft2(x12);

  //   Eigen::MatrixXd x33(3, 3);
  //   x33 << 2, -1.4, 1, -9, 2, 3.9, 13, 1.3, -2.2;
  //   expect_rfft2(x33);
}

// void expect_inv_rfft2(const Eigen::MatrixXcd& x) {
//   for (int n = 0; n < x.cols(); ++n) {
//     for (int m = 0; m < x.rows(); ++m) {
//       auto g = [m, n](const auto& x) {
//         using stan::math::inv_rfft2;
//         return inv_rfft2(x)(m, n);
//       };
//       stan::test::expect_ad(g, x);
//     }
//   }
// }

// TEST(mathMixFun, inv_fft2) {
//   using cmat_t = Eigen::MatrixXcd;

//   cmat_t x00(0, 0);
//   expect_inv_fft2(x00);

//   cmat_t x11(1, 1);
//   x11(0, 0) = {1, 2};
//   expect_inv_fft2(x11);

//   cmat_t x21(2, 1);
//   x21(0, 0) = {1, 2};
//   x21(1, 0) = {-1.3, 2.9};
//   expect_inv_fft2(x21);

//   cmat_t x22(2, 2);
//   x22(0, 0) = {3, 9};
//   x22(0, 1) = {-13.2, 8.345};
//   x22(1, 0) = {-4.23, 7.566};
//   x22(1, 1) = {1, -12.9};
//   expect_inv_fft2(x22);

//   cmat_t x33(3, 3);
//   x33(0, 0) = {3, 9};
//   x33(0, 1) = {-13.2, 8.345};
//   x33(0, 2) = {4.333, -1.9};
//   x33(1, 0) = {-4.23, 7.566};
//   x33(1, 1) = {1, -12.9};
//   x33(1, 2) = {-1.01, -4.01};
//   x33(2, 0) = {3.87, 5.89};
//   x33(2, 1) = {-2.875, -2.999};
//   x33(2, 2) = {12.98, 14.5555};
//   expect_inv_fft2(x33);
// }
