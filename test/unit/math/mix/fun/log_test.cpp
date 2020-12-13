#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(mathMixMatFun, log) {
  auto f = [](const auto& x1) {
    using stan::math::log;
    return log(x1);
  };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -0.2, 1e-3, 1, 1.3, 3, 3.7, 10, 10.2,
                                      1e6);

  // non-zero real and imaginary components
  for (auto re : std::vector<double>{-2.7, 1, 2.3}) {
    for (auto im : std::vector<double>{-1.5, 1.2}) {
      stan::test::expect_ad(f, std::complex<double>{re, im});
    }
  }
  // zero tests which lead to finite vals
  stan::test::expect_ad(f, std::complex<double>{0, 2.1});
  stan::test::expect_ad(f, std::complex<double>{0, -2.1});
  stan::test::expect_ad(f, std::complex<double>{2.1, 0});
  stan::test::expect_ad(f, std::complex<double>{-0.0, 2.1});
  stan::test::expect_ad(f, std::complex<double>{-0.0, -2.1});
  stan::test::expect_ad(f, std::complex<double>{2.1, -0.0});
  // (negative real and zero imaginary illegal)

  Eigen::VectorXd x1(3);
  x1 << 0.1, 2.0, 3.0;
  Eigen::RowVectorXd x2(3);
  x2 << 0.1, 2.0, 3.0;
  Eigen::MatrixXd x3(2, 3);
  x3 << 0.1, 2.0, 3.0, 4.0, 5.0, 6.0;
  stan::test::expect_ad_matvar(f, x1);
  stan::test::expect_ad_matvar(f, x2);
  stan::test::expect_ad_matvar(f, x3);

  Eigen::VectorXd x4(0);
  Eigen::RowVectorXd x5(0);
  Eigen::MatrixXd x6(0, 0);

  stan::test::expect_ad_matvar(f, x4);
  stan::test::expect_ad_matvar(f, x5);
  stan::test::expect_ad_matvar(f, x6);
}
