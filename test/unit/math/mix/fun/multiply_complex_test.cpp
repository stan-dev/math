#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/fun/multiply_util.hpp>

TEST(mathMix, complexMultiply) {
  auto f
      = [](const auto& x, const auto& y) { return stan::math::multiply(x, y); };
  double d_scalar = 1.0;
  std::complex<double> cd_scalar = 1.0;
  Eigen::Matrix<double, -1, -1> d_mat(2, 2);
  d_mat << 1, 2, 3, 4;
  Eigen::Matrix<std::complex<double>, -1, -1> cd_mat(2, 2);
  cd_mat << 1, 2, 3, 4;

  stan::test::expect_ad(f, d_scalar, cd_mat);
  stan::test::expect_ad(f, cd_scalar, d_mat);
  stan::test::expect_ad(f, cd_scalar, cd_mat);
  stan::test::expect_ad(f, cd_mat, d_scalar);
  stan::test::expect_ad(f, d_mat, cd_scalar);
  stan::test::expect_ad(f, cd_mat, cd_scalar);
  stan::test::expect_ad(f, d_mat, cd_mat);
  stan::test::expect_ad(f, cd_mat, d_mat);
  stan::test::expect_ad(f, cd_mat, cd_mat);

  Eigen::Matrix<double, -1, 1> d_vec(2);
  d_vec << 1, 2;
  Eigen::Matrix<std::complex<double>, -1, 1> cd_vec(2);
  cd_vec << 1, 2;

  stan::test::expect_ad(f, d_mat, cd_vec);
  stan::test::expect_ad(f, cd_mat, d_vec);
  stan::test::expect_ad(f, cd_mat, cd_vec);

  Eigen::Matrix<double, 1, -1> d_rowvec(2);
  d_rowvec << 1, 2;
  Eigen::Matrix<std::complex<double>, 1, -1> cd_rowvec(2);
  cd_rowvec << 1, 2;

  stan::test::expect_ad(f, d_vec, cd_rowvec);
  stan::test::expect_ad(f, cd_vec, d_rowvec);
  stan::test::expect_ad(f, cd_vec, cd_rowvec);
  stan::test::expect_ad(f, d_rowvec, cd_vec);
  stan::test::expect_ad(f, cd_rowvec, d_vec);
  stan::test::expect_ad(f, cd_rowvec, cd_vec);
  stan::test::expect_ad(f, d_rowvec, cd_mat);
  stan::test::expect_ad(f, cd_rowvec, d_mat);
  stan::test::expect_ad(f, cd_rowvec, cd_mat);
  stan::test::expect_ad(f, d_mat, cd_vec);
  stan::test::expect_ad(f, cd_mat, d_vec);
  stan::test::expect_ad(f, cd_mat, cd_vec);
}
