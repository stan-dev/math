#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

TEST(AgradRevMatrix, eigenval_sum) {
  using stan::math::sum;
  using stan::math::matrix_v;
  using stan::math::vector_v;
  using stan::math::eigenvalues_sym;

  matrix_v a(3, 3);
  a <<
    1.0, 2.0, 3.0,
    2.0, 5.0, 7.9,
    3.0, 7.9, 1.08;
  AVEC x = createAVEC(a(0, 0), a(1, 1), a(2, 2), a(1, 2));
  x.push_back(a(0, 1));
  x.push_back(a(2, 0));

  // grad sum eig = I
  vector_v a_eigenvalues = eigenvalues_sym(a);

  AVAR sum_a_eigenvalues = sum(a_eigenvalues);

  VEC g = cgradvec(sum_a_eigenvalues, x);

  EXPECT_NEAR(1.0, g[0], 1.0E-11);
  EXPECT_NEAR(1.0, g[1], 1.0E-11);
  EXPECT_NEAR(1.0, g[2], 1.0E-11);

  EXPECT_NEAR(0.0, g[3], 1.0E-10);
  EXPECT_NEAR(0.0, g[4], 1.0E-10);
  EXPECT_NEAR(0.0, g[5], 1.0E-10);
}
TEST(AgradRevMatrix, check_varis_on_stack) {
  using stan::math::matrix_v;
  using stan::math::eigenvalues_sym;

  matrix_v a(3, 3);
  a <<
    1.0, 2.0, 3.0,
    2.0, 5.0, 7.9,
    3.0, 7.9, 1.08;
  test::check_varis_on_stack(stan::math::eigenvalues_sym(a));
}
