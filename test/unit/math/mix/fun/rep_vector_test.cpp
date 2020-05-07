#include <test/unit/math/test_ad.hpp>
namespace rep_vector_test { 
auto f(int n) {
  return [=](const auto& y) { return stan::math::rep_vector(y, n); };
}
}

TEST(MathMixMatFun, repVector) {
  double y = 3;
  stan::test::expect_ad(rep_vector_test::f(0), y);
  stan::test::expect_ad(rep_vector_test::f(1), y);
  stan::test::expect_ad(rep_vector_test::f(4), y);
}
