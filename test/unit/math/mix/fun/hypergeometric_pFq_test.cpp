#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mixScalFun, grad_2F2_ffv) {
  using stan::math::fvar;
  using stan::math::hypergeometric_pFq;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_ffv;

  vector_ffv ffv_a(2);
  vector_d d_a(2);
  ffv_a.val().val() << 4, 2;
  d_a << 4, 2;
  ffv_a.val().d() << 1, 1;

  vector_ffv ffv_b(2);
  vector_d d_b(2);
  ffv_b.val().val() << 6, 3;
  d_b << 6, 3;
  ffv_b.val().d() << 1, 1;

  fvar<fvar<var>> ffv_z;
  ffv_z.val_.val_ = 4;
  ffv_z.val_.d_ = 1;
  double d_z = 4;

  double a_adj = 3.924636646666071 + 6.897245961898751;
  double b_adj = -2.775051002566842 - 4.980095849781222;
  double z_adj = 4.916522138006060;

  // fvar, fvar, fvar
  EXPECT_FLOAT_EQ(hypergeometric_pFq(ffv_a, ffv_b, ffv_z).val_.d_.val(),
                  a_adj + b_adj + z_adj);
  // fvar, fvar, double
  EXPECT_FLOAT_EQ(hypergeometric_pFq(ffv_a, ffv_b, d_z).val_.d_.val(),
                  a_adj + b_adj);
  // fvar, double, double
  EXPECT_FLOAT_EQ(hypergeometric_pFq(ffv_a, d_b, d_z).val_.d_.val(), a_adj);
  // fvar, double, fvar
  EXPECT_FLOAT_EQ(hypergeometric_pFq(ffv_a, d_b, ffv_z).val_.d_.val(),
                  a_adj + z_adj);
  // double, fvar, fvar
  EXPECT_FLOAT_EQ(hypergeometric_pFq(d_a, ffv_b, ffv_z).val_.d_.val(),
                  b_adj + z_adj);
  // double, fvar, double
  EXPECT_FLOAT_EQ(hypergeometric_pFq(d_a, ffv_b, d_z).val_.d_.val(), b_adj);
  // double, double, fvar
  EXPECT_FLOAT_EQ(hypergeometric_pFq(d_a, d_b, ffv_z).val_.d_.val(), z_adj);
}
