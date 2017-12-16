#include <stan/math/rev/mat.hpp>
#include <stan/math/prim/mat/fun/matrix_exp_2x2.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/expect_matrix_eq.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>

TEST(MathMatrix, matrix_exp_2x2) {
  // example from Moler & Van Loan, 2003
  for (size_t k = 0; k < 2; k++) {
    for (size_t l = 0; l < 2; l++) {
      AVAR a = -1.0, b = -17.0;

      stan::math::matrix_v m1(2, 2), m2(2, 2), m1_exp;
      m1 << -2 * a + 3 * b, 1.5 * a - 1.5 * b, -4 * a + 4 * b, 3 * a - 2 * b;
      m2 << -.735759, .551819, -1.471518, 1.103638;

      m1_exp = stan::math::matrix_exp_2x2(m1);
      expect_matrix_eq(m2, m1_exp);

      stan::math::matrix_v dm1_exp_da(2, 2), dm1_exp_db(2, 2);
      dm1_exp_da << -2 * exp(a), 1.5 * exp(a), -4 * exp(a), 3 * exp(a);
      dm1_exp_db << 3 * exp(b), -1.5 * exp(b), 4 * exp(b), -2 * exp(b);

      AVEC x = createAVEC(a, b);
      VEC g;
      m1_exp(k, l).grad(x, g);
      EXPECT_FLOAT_EQ(dm1_exp_da(k, l).val(), g[0]);
      EXPECT_FLOAT_EQ(dm1_exp_db(k, l).val(), g[1]);
    }
  }
}
