#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(primScalFun, grad_2F2_ffv) {
  using stan::math::vector_ffv;
  using stan::math::vector_d;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::hypergeometric_pFq;

  vector_ffv ffv_p(2);
  vector_d d_p(2);
  ffv_p.val() << 4, 2;
  d_p << 4, 2;
  ffv_p.val().d() << 1, 1;

  vector_ffv ffv_q(2);
  vector_d d_q(2);
  ffv_q.val() << 6, 3;
  d_q << 6, 3;
  ffv_q.val().d() << 1, 1;

  fvar<fvar<var>> ffv_z;
  ffv_z.val_.val_ = 4;
  ffv_z.val_.d_ = 1;
  double d_z = 4;

  double p_adj = 3.924636646666071 + 6.897245961898751;
  double q_adj = -2.775051002566842 - 4.980095849781222;
  double z_adj = 4.916522138006060;

  // fvar, fvar, fvar
  EXPECT_FLOAT_EQ(hypergeometric_pFq(ffv_p, ffv_q, ffv_z).val_.d_.val(),
                  p_adj + q_adj + z_adj);
  // fvar, fvar, double
  EXPECT_FLOAT_EQ(hypergeometric_pFq(ffv_p, ffv_q, d_z).val_.d_.val(),
                  p_adj + q_adj);
  // fvar, double, double
  EXPECT_FLOAT_EQ(hypergeometric_pFq(ffv_p, d_q, d_z).val_.d_.val(), p_adj);
  // fvar, double, fvar
  EXPECT_FLOAT_EQ(hypergeometric_pFq(ffv_p, d_q, ffv_z).val_.d_.val(),
                  p_adj + z_adj);
  // double, fvar, fvar
  EXPECT_FLOAT_EQ(hypergeometric_pFq(d_p, ffv_q, ffv_z).val_.d_.val(),
                  q_adj + z_adj);
  // double, fvar, double
  EXPECT_FLOAT_EQ(hypergeometric_pFq(d_p, ffv_q, d_z).val_.d_.val(), q_adj);
  // double, double, fvar
  EXPECT_FLOAT_EQ(hypergeometric_pFq(d_p, d_q, ffv_z).val_.d_.val(), z_adj);
}
