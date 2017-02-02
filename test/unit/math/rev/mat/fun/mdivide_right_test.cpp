#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

TEST(AgradRevMatrix,mdivide_right_val) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::mdivide_right;

  matrix_v Av(2,2);
  matrix_d Ad(2,2);
  matrix_v I;

  Av << 2.0, 3.0, 
    5.0, 7.0;
  Ad << 2.0, 3.0, 
    5.0, 7.0;

  I = mdivide_right(Av,Av);
  EXPECT_NEAR(1.0,I(0,0).val(),1.0E-12);
  EXPECT_NEAR(0.0,I(0,1).val(),1.0E-12);
  EXPECT_NEAR(0.0,I(1,0).val(),1.0E-12);
  EXPECT_NEAR(1.0,I(1,1).val(),1.0e-12);

  I = mdivide_right(Av,Ad);
  EXPECT_NEAR(1.0,I(0,0).val(),1.0E-12);
  EXPECT_NEAR(0.0,I(0,1).val(),1.0E-12);
  EXPECT_NEAR(0.0,I(1,0).val(),1.0E-12);
  EXPECT_NEAR(1.0,I(1,1).val(),1.0e-12);

  I = mdivide_right(Ad,Av);
  EXPECT_NEAR(1.0,I(0,0).val(),1.0E-12);
  EXPECT_NEAR(0.0,I(0,1).val(),1.0E-12);
  EXPECT_NEAR(0.0,I(1,0).val(),1.0E-12);
  EXPECT_NEAR(1.0,I(1,1).val(),1.0e-12);
}

TEST(AgradRevMatrix, check_varis_on_stack) {
  using stan::math::value_of;
  stan::math::matrix_v m(2,2);
  m << 2.0, 3.0, 
    5.0, 7.0;

  test::check_varis_on_stack(stan::math::mdivide_right(m, m));
  test::check_varis_on_stack(stan::math::mdivide_right(m, value_of(m)));
  test::check_varis_on_stack(stan::math::mdivide_right(value_of(m), m));
}
