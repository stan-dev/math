#include <test/unit/math/test_ad.hpp>
#include <cmath>
#include <complex>
#include <vector>
#include <type_traits>

TEST(mixFun, absBasics) {
  using stan::math::abs;
  int a = abs(1);

  double b = abs(-2.3);

  Eigen::Matrix<double, -1, -1> x(2, 3);
  x << 1, 2, 3, 4, 5, 6;
  Eigen::Matrix<double, -1, -1> y = abs(x);

  std::vector<int> u{1, 2, 3, 4};
  std::vector<int> v = abs(u);
}

TEST(mixFun, abs) {
  auto f = [](const auto& x) { return stan::math::abs(x); };
  stan::test::expect_common_nonzero_unary(f);
  // 0 (no derivative at 0)
  stan::test::expect_value(f, 0);
  stan::test::expect_value(f, 0.0);

  stan::test::expect_ad(f, -3);
  stan::test::expect_ad(f, -2);
  stan::test::expect_ad(f, 2);

  stan::test::expect_ad(f, -17.3);
  stan::test::expect_ad(f, -0.68);
  stan::test::expect_ad(f, 0.68);
  stan::test::expect_ad(f, 2.0);
  stan::test::expect_ad(f, 4.0);

  // complex tests
  for (double re : std::vector<double>{-4, -2.5, -1.5, -0.3, 1.3, 2.1, 3.9}) {
    for (double im : std::vector<double>{-4, -2.5, -1.5, -0.3, 1.3, 2.1, 3.9}) {
      stan::test::expect_ad(f, std::complex<double>(re, im));
    }
  }

  // vector<double>
  using svd_t = std::vector<double>;
  stan::test::expect_ad(f, svd_t{});
  stan::test::expect_ad(f, svd_t{1.0});
  stan::test::expect_ad(f, svd_t{1.9, -2.3});

  // vector<vector<double>>
  using svvd_t = std::vector<svd_t>;
  stan::test::expect_ad(f, svvd_t{});
  stan::test::expect_ad(f, svvd_t{svd_t{}});
  stan::test::expect_ad(f, svvd_t{svd_t{1.9, 4.8}});
  stan::test::expect_ad(f, svvd_t{svd_t{1.9}, svd_t{-13.987}});
  stan::test::expect_ad(f, svvd_t{svd_t{1.9, -2.7}, svd_t{-13.987, 8.8}});

  // vector<complex<double>>
  using c_t = std::complex<double>;
  using svc_t = std::vector<c_t>;
  stan::test::expect_ad(f, svc_t{});
  stan::test::expect_ad(f, svc_t{c_t{1.0, -1.9}});
  stan::test::expect_ad(f, svc_t{c_t{1.0, -1.9}, c_t{-9.3, -128.987654}});

  // vector<vector<complex<double>>>
  using svvc_t = std::vector<svc_t>;
  stan::test::expect_ad(f, svvc_t{});
  stan::test::expect_ad(f, svvc_t{{}});
  stan::test::expect_ad(f, svvc_t{svc_t{c_t{1.2, -2.3}, c_t{-32.8, 1}}});
  stan::test::expect_ad(f, svvc_t{svc_t{c_t{1.2, -2.3}, c_t{-32.8, 1}},
                                  svc_t{c_t{9.3, 9.4}, c_t{182, -95}}});

  // VectorXd
  using v_t = Eigen::VectorXd;
  v_t a0(0);
  stan::test::expect_ad(f, a0);
  stan::test::expect_ad_matvar(f, a0);
  v_t a1(1);
  a1 << 1.9;
  stan::test::expect_ad(f, a1);
  stan::test::expect_ad_matvar(f, a1);
  v_t a2(2);
  a2 << 1.9, -2.3;
  stan::test::expect_ad(f, a2);
  stan::test::expect_ad_matvar(f, a2);

  // RowVectorXd
  using rv_t = Eigen::RowVectorXd;
  rv_t b0(0);
  stan::test::expect_ad(f, b0);
  stan::test::expect_ad_matvar(f, b0);
  rv_t b1(1);
  b1 << 1.9;
  stan::test::expect_ad(f, b1);
  stan::test::expect_ad_matvar(f, b1);
  rv_t b2(2);
  b2 << 1.9, -2.3;
  stan::test::expect_ad(f, b2);
  stan::test::expect_ad_matvar(f, b2);

  // MatrixXd
  using m_t = Eigen::MatrixXd;
  m_t c0(0, 0);
  stan::test::expect_ad(f, c0);
  stan::test::expect_ad_matvar(f, c0);
  m_t c0i(0, 2);
  stan::test::expect_ad(f, c0i);
  stan::test::expect_ad_matvar(f, c0i);
  m_t c0ii(2, 0);
  stan::test::expect_ad(f, c0ii);
  stan::test::expect_ad_matvar(f, c0ii);
  m_t c2(2, 1);
  c2 << 1.3, -2.9;
  stan::test::expect_ad(f, c2);
  stan::test::expect_ad_matvar(f, c2);
  m_t c6(3, 2);
  c6 << 1.3, 2.9, -13.456, 1.898, -0.01, 1.87e21;
  stan::test::expect_ad(f, c6);
  stan::test::expect_ad_matvar(f, c6);

  // vector<VectorXd>
  using av_t = std::vector<Eigen::VectorXd>;
  av_t d0;
  stan::test::expect_ad(f, d0);
  stan::test::expect_ad_matvar(f, d0);
  av_t d1{a0};
  stan::test::expect_ad(f, d1);
  stan::test::expect_ad_matvar(f, d1);
  av_t d2{a1, a2};
  stan::test::expect_ad(f, d2);
  stan::test::expect_ad_matvar(f, d2);

  // vector<RowVectorXd>
  using arv_t = std::vector<Eigen::RowVectorXd>;
  arv_t e0;
  stan::test::expect_ad(f, e0);
  stan::test::expect_ad_matvar(f, e0);
  arv_t e1{b0};
  stan::test::expect_ad(f, e1);
  stan::test::expect_ad_matvar(f, e1);
  arv_t e2{b1, b2};
  stan::test::expect_ad(f, e2);
  stan::test::expect_ad_matvar(f, e2);

  // vector<MatrixXd>
  using am_t = std::vector<Eigen::MatrixXd>;
  am_t g0;
  stan::test::expect_ad(f, g0);
  stan::test::expect_ad_matvar(f, g0);
  am_t g1{c0};
  stan::test::expect_ad(f, g1);
  stan::test::expect_ad_matvar(f, g1);
  am_t g2{c2, c6};
  stan::test::expect_ad(f, g2);
  stan::test::expect_ad_matvar(f, g2);

  // VectorXcd
  using vc_t = Eigen::VectorXcd;
  vc_t h0(0);
  stan::test::expect_ad(f, h0);
  vc_t h1(1);
  h1 << c_t{1.9, -1.8};
  stan::test::expect_ad(f, h1);
  vc_t h2(2);
  h2 << c_t{1.9, -1.8}, c_t{-128.7, 1.3};
  stan::test::expect_ad(f, h2);

  // RowVectorXcd
  using rvc_t = Eigen::RowVectorXcd;
  rvc_t j0(0);
  stan::test::expect_ad(f, j0);
  rvc_t j1(1);
  j1 << c_t{1.9, -1.8};
  stan::test::expect_ad(f, j1);
  rvc_t j2(2);
  j2 << c_t{1.9, -1.8}, c_t{-128.7, 1.3};
  stan::test::expect_ad(f, j2);

  // MatrixXcd
  using mc_t = Eigen::MatrixXcd;
  mc_t k0(0, 0);
  stan::test::expect_ad(f, k0);
  mc_t k2(1, 2);
  k2 << c_t{1.9, -1.8}, c_t{128.735, 128.734};
  stan::test::expect_ad(f, k2);
  mc_t k6(3, 2);
  k6 << c_t{1.9, -1.8}, c_t{-128.7, 1.3}, c_t{1, 2}, c_t{0.3, -0.5},
      c_t{-13, 125.7}, c_t{-12.5, -10.5};
  stan::test::expect_ad(f, k6);

  // vector<VectorXcd>
  using avc_t = std::vector<vc_t>;
  avc_t m0;
  stan::test::expect_ad(f, m0);
  avc_t m1{h1};
  stan::test::expect_ad(f, m1);
  avc_t m2{h1, h2};
  stan::test::expect_ad(f, m2);

  // vector<RowVectorXcd>
  using arvc_t = std::vector<rvc_t>;
  arvc_t p0(0);
  stan::test::expect_ad(f, p0);
  arvc_t p1{j1};
  stan::test::expect_ad(f, p1);
  arvc_t p2{j1, j2};
  stan::test::expect_ad(f, p2);

  // vector<MatrixXcd>
  using amc_t = std::vector<mc_t>;
  amc_t q0;
  stan::test::expect_ad(f, q0);
  amc_t q1{k2};
  stan::test::expect_ad(f, q1);
  amc_t q2{k2, k6};
  stan::test::expect_ad(f, q2);
}
TEST(mixFun, absReturnType) {
  // validate return types not overpromoted to complex by assignability
  std::complex<stan::math::var> a = 3;
  stan::math::var b = abs(a);

  std::complex<stan::math::fvar<double>> c = 3;
  stan::math::fvar<double> d = abs(c);
  SUCCEED();
}

TEST(mathMixMatFun, abs_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_nonzero_args;
  auto f = [](const auto& x1) {
    using stan::math::abs;
    return abs(x1);
  };
  std::vector<double> com_args = common_nonzero_args();
  std::vector<double> args{-3, 2, -0.68, 1};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}
