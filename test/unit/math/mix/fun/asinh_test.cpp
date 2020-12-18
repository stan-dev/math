#include <test/unit/math/test_ad.hpp>
#include <complex>
#include <vector>

TEST(mathMixFun, asinh) {
  auto f = [](const auto& x1) {
    using stan::math::asinh;
    return asinh(x1);
  };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -1.2, -0.2, 0.5, 2, -1.2);
  stan::test::expect_ad(f, std::complex<double>{0.2});
  // avoid pole at real zero that can't be autodiffed
  for (double re : std::vector<double>{-0.2, 0.3}) {
    for (double im : std::vector<double>{-0.3, 0, 0.2}) {
      stan::test::expect_ad(f, std::complex<double>{re, im});
    }
  }
}

TEST(mathMixMatFun, asinh_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_args;
  auto f = [](const auto& x1) {
    using stan::math::asinh;
    return asinh(x1);
  };
  std::vector<double> com_args = common_args();
  std::vector<double> args{-2.6, -1.2, -0.2, 0.5, 2, -1.2};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}
