#include <test/unit/math/test_ad.hpp>
#include <array>
#include <limits>
#include <vector>

TEST(mathMixMatFun, sech) {
  auto f = [](auto&& x1) { return stan::math::sech(x1); };
  stan::test::expect_common_unary_vectorized<stan::test::ScalarSupport::Real>(
      f);
  stan::test::expect_unary_vectorized<stan::test::ScalarSupport::Real>(
      f, -30.0, -2.6, -1.2, -0.2, 0.0, 0.5, 1.3, 30.0);
}

TEST(mathMixMatFun, sech_varmat) {
  using stan::math::vec_concat;
  using stan::test::expect_ad_vector_matvar;
  using stan::test::internal::common_args;
  auto f = [](auto&& x1) { return stan::math::sech(x1); };
  std::vector<double> com_args = common_args();
  std::vector<double> args{-30.0, -2.6, -1.2, -0.2, 0.5, 1.3, 30.0};
  auto all_args = vec_concat(com_args, args);
  Eigen::VectorXd A(all_args.size());
  for (int i = 0; i < all_args.size(); ++i) {
    A(i) = all_args[i];
  }
  expect_ad_vector_matvar(f, A);
}

// The var_value Cholesky factor builds exp(-acc) as a product of sech; it
// must stay accurate through its underflow at 800.
TEST(mathMixMatFun, sech_reference) {
  const double eps = std::numeric_limits<double>::epsilon();
  // x, sech(x)
  std::vector<std::array<double, 2>> cases = {{0.0, 1.0},
                                              {1e-10, 1.0},
                                              {0.5, 0.88681888397007391},
                                              {0.999, 0.64854787974731328},
                                              {1.0, 0.6480542736638854},
                                              {1.001, 0.64756077130222642},
                                              {-3.0, 0.099327927419433208},
                                              {30.0, 1.8715245937680349e-13},
                                              {700.0, 1.9719353087519542e-304},
                                              {800.0, 0.0}};
  for (auto&& c : cases) {
    EXPECT_NEAR(c[1], stan::math::sech(c[0]), 4 * eps * c[1]) << "x = " << c[0];
  }
}
