#include <test/unit/math/test_ad.hpp>
namespace rep_row_vector_test { 
auto f(int n) {
  return [=](const auto& y) { return stan::math::rep_row_vector(y, n); };
}
}

TEST(MathMixMatFun, repRowVector) {
  double y = 3;
  stan::test::expect_ad(rep_row_vector_test::f(0), y);
  stan::test::expect_ad(rep_row_vector_test::f(1), y);
  stan::test::expect_ad(rep_row_vector_test::f(4), y);
}
