#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, cos) {
  auto f = [](const auto& x) {
    using stan::math::cos;
    return cos(x);
  };
  stan::test::expect_common_nonzero_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2, -0.2, -0.5, 0, 1.5, 3, 5,
                                      5.3);
  stan::test::expect_complex_common(f);
}

TEST(mathMixMatFun, cos_varmat) {
  auto f = [](const auto& x1) {
    using stan::math::cos;
    return cos(x1);
  };
  auto com_args = stan::test::internal::common_nonzero_args();
  std::vector<double> extra_args{-2.6, -2, -0.2, -0.5, 0, 1.5, 3, 5, 5.3};
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
