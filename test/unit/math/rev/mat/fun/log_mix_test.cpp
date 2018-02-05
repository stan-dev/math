#include <stan/math/rev/mat.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <gtest/gtest.h>
#include <vector>

using stan::math::log_mix;
using stan::math::vector_d;
using stan::math::vector_v;
using stan::math::row_vector_v;
using stan::math::var;

template <typename T_a, typename T_b>
void val_rev_test(T_a a, T_b b) {
  a[0] = 0.112;
  a[1] = 0.214;
  a[2] = 0.305;
  a[3] = 0.369;

  b[0] = -5.983;
  b[1] = -11.215;
  b[2] = -6.836;
  b[3] = -4.538;

  AVAR out = log_mix(a, b);

  out.grad();

  EXPECT_FLOAT_EQ(out.val(), -5.390580825249);

  vector_d prob_deriv(4);
  prob_deriv << 0.55298789093, 0.00295451971, 0.23564727846, 2.3456928701;
  vector_d dens_deriv(4);
  dens_deriv << 0.06193464378, 0.00063226722, 0.07187241993, 0.86556066907;


  for (int i = 0; i < 4; ++i) {
    EXPECT_FLOAT_EQ(a[i].adj(), prob_deriv[i]);
    EXPECT_FLOAT_EQ(b[i].adj(), dens_deriv[i]);
  }
};

template <typename T_a, typename T_b>
void val_rev_vec_test(T_a a, T_b b1) {
  a[0] = 0.112;
  a[1] = 0.214;
  a[2] = 0.305;
  a[3] = 0.369;

  T_b b2(4);
  T_b b3(4);

  b1[0] = -5.983;
  b1[1] = -11.215;
  b1[2] = -6.836;
  b1[3] = -4.538;

  b2[0] = -10.365;
  b2[1] = -12.443;
  b2[2] = -15.091;
  b2[3] = -19.115;

  b3[0] = -4.117;
  b3[1] = -8.132;
  b3[2] = -7.931;
  b3[3] = -12.115;

  std::vector<T_b> c{b1, b2, b3};

  AVAR out = log_mix(a, c);

  out.grad();

  EXPECT_FLOAT_EQ(out.val(), -23.9255869110932);

  vector_d prob_deriv(4);
  prob_deriv << 15.7666988503552, 1.03434231079226, 0.478019548364374, 2.34955152303399;
  vector_d dens1_deriv(4);
  dens1_deriv << 0.061934643784246, 0.000632267218743, 0.071872419930696, 0.865560669066315;
  vector_d dens2_deriv(4);
  dens2_deriv << 0.791240264915128, 0.189251877628038, 0.019094771903811, 0.000413085553023;
  vector_d dens3_deriv(4);
  dens3_deriv << 0.91269536254040, 0.031465109662762, 0.054828770416628, 0.001010757380203;

  for (int i = 0; i < 4; ++i) {
    EXPECT_FLOAT_EQ(a[i].adj(), prob_deriv[i]);
    EXPECT_FLOAT_EQ(c[0][i].adj(), dens1_deriv[i]);
    EXPECT_FLOAT_EQ(c[1][i].adj(), dens2_deriv[i]);
    EXPECT_FLOAT_EQ(c[2][i].adj(), dens3_deriv[i]);
  }
};

TEST(AgradRevMatrix, log_mix_vals) {
  vector_v vecv_prob(4), vecv_dens(4);
  row_vector_v row_vecv_prob(4), row_vecv_dens(4);
  std::vector<var> stvec_prob(4), stvec_dens(4);

  val_rev_test(vecv_prob, vecv_dens);
  val_rev_test(vecv_prob, row_vecv_dens);
  val_rev_test(vecv_prob, stvec_dens);
  val_rev_vec_test(vecv_prob, vecv_dens);
  val_rev_vec_test(vecv_prob, row_vecv_dens);

  val_rev_test(row_vecv_prob, vecv_dens);
  val_rev_test(row_vecv_prob, row_vecv_dens);
  val_rev_test(row_vecv_prob, stvec_dens);
  val_rev_vec_test(row_vecv_prob, vecv_dens);
  val_rev_vec_test(row_vecv_prob, row_vecv_dens);

  val_rev_test(stvec_prob, vecv_dens);
  val_rev_test(stvec_prob, row_vecv_dens);
  val_rev_test(stvec_prob, stvec_dens);
  val_rev_vec_test(stvec_prob, vecv_dens);
  val_rev_vec_test(stvec_prob, row_vecv_dens);
}
