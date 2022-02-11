#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(mathMixMatFun, atanh) {
  auto f = [](const auto& x1) {
    using stan::math::atanh;
    return atanh(x1);
  };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -0.9, 0.5);
  for (double re : std::vector<double>{-0.2, 0, 0.3}) {
    for (double im : std::vector<double>{-0.3, 0, 0.2}) {
      stan::test::expect_ad(f, std::complex<double>{re, im});
    }
  }
}

TEST(mathMixMatFun, atanh_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_args;
  auto f = [](const auto& x1) {
    using stan::math::atanh;
    return atanh(x1);
  };
  std::vector<double> com_args = common_args();
  std::vector<double> args{-0.9, 0.5};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}
