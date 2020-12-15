#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, cosh) {
  auto f = [](const auto& x1) {
    using stan::math::cosh;
    return cosh(x1);
  };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2, -1.2, -0.2, 0.5, 1, 1.3,
                                      1.5);
  stan::test::expect_complex_common(f);
}

TEST(mathMixMatFun, cosh_varmat) {
  auto f = [](const auto& x1) {
    using stan::math::cosh;
    return cosh(x1);
  };
  auto com_args = stan::test::internal::common_args();
  std::vector<double> extra_args{-2.6, -2, -1.2, -0.2, 0.5, 1, 1.3, 1.5};
  Eigen::VectorXd A(com_args.size() + extra_args.size());
  int i = 0;
  for (double x : com_args) {
    A(i) = x;
    ++i;
  }
  for (double x : extra_args) {
    A(i) = x;
    ++i;
  }
  stan::test::expect_ad_matvar(f, A);
  std::vector<Eigen::VectorXd> A_vec;
  A_vec.push_back(A);
  A_vec.push_back(A);
  A_vec.push_back(A);
  stan::test::expect_ad_matvar(f, A_vec);
}
