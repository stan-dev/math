#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(mathMixMatFun, acosh) {
  auto f = [](const auto& x1) {
    using stan::math::acosh;
    return acosh(x1);
  };
  for (double x : stan::test::internal::common_args())
    stan::test::expect_unary_vectorized(x);
  stan::test::expect_unary_vectorized(f, 1.5, 3.2, 5, 10, 12.9);
  // avoid pole at complex zero that can't be autodiffed
  for (double re : std::vector<double>{-0.2, 0, 0.3}) {
    for (double im : std::vector<double>{-0.3, 0.2}) {
      stan::test::expect_ad(f, std::complex<double>{re, im});
    }
  }
}

TEST(mathMixMatFun, acosh_varmat) {
  auto f = [](const auto& x1) {
    using stan::math::acosh;
    return acosh(x1);
  };
  auto com_args = stan::test::internal::common_args();
  std::vector<double> extra_args{1.5, 3.2, 5, 10, 12.9};
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
  std::vector<Eigen::MatrixXd> A_vec;
  A_vec.push_back(A);
  A_vec.push_back(A);
  A_vec.push_back(A);
  stan::test::expect_ad_matvar(f, A_vec);

}
