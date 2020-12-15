#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, fabs) {
  auto f = [](const auto& x1) { return stan::math::fabs(x1); };
  stan::test::expect_common_nonzero_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2, -0.5, 1.5, 2.0, 3);
}

TEST(mathMixMatFun, fabs_varmat) {
  auto f = [](const auto& x1) {
    using stan::math::fabs;
    return fabs(x1);
  };
  auto com_args = stan::test::internal::common_args();
  std::vector<double> extra_args{-2.6, -2, -0.5, 1.5, 2.0, 3};
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
