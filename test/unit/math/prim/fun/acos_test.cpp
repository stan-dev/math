#include <stan/math/prim.hpp>
#include <test/unit/math/test_ad.hpp>

TEST(primScalFun, acos) {
  stan::test::expect_common_prim([](auto x) { return std::acos(x); },
                                 [](auto x) { return stan::math::acos(x); });
}

TEST(MathFunctions, acos_works_with_other_functions) {
  Eigen::VectorXd a(5);
  a << 0.1, 0.2, 0.3, 0.4, 0.5;
  Eigen::RowVectorXd b(5);
  b << 0.1, 0.2, 0.3, 0.4, 0.5;
  stan::math::multiply(a, stan::math::acos(b));
}
