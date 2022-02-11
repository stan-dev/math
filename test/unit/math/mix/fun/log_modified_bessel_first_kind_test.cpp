#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(mathMixScalFun, logModifiedBesselFirstKind) {
  auto f = [](const auto& x1, const auto& x2) {
    return stan::math::log_modified_bessel_first_kind(x1, x2);
  };

  std::vector<double> xs{-1.3, 0.5, stan::math::positive_infinity(),
                         stan::math::negative_infinity(),
                         stan::math::not_a_number()};

  for (double x1 : xs)
    for (double x2 : xs)
      stan::test::expect_ad(f, x1, x2);
}

TEST(mathMixScalFun, logModifiedBesselFirstKind_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::log_modified_bessel_first_kind;
    return log_modified_bessel_first_kind(x1, x2);
  };

  Eigen::VectorXd in1(2);
  in1 << -1.3, 0.5;
  Eigen::VectorXd in2(2);
  in2 << 0.5, 3.4;
  stan::test::expect_ad_vectorized_binary(f, in1, in2);
}
