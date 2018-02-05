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
  a[0] = 0.15;
  a[1] = 0.20;
  a[2] = 0.40;
  a[3] = 0.25;

  b[0] = -2.15;
  b[1] = -3.89;
  b[2] = -2.18;
  b[3] = -8.82;

  AVAR out = log_mix(a, b);

  out.grad();

  EXPECT_FLOAT_EQ(out.val(), -2.70582405);

  vector_d prob_deriv(4);
  prob_deriv << 1.7433770185, 0.3059982327, 1.6918524409, 0.0022112972;
  vector_d dens_deriv(4);
  dens_deriv << 0.2615065527, 0.0611996465, 0.6767409763, 0.0005528243;


  for (int i = 0; i < 4; ++i) {
    EXPECT_FLOAT_EQ(a[i].adj(), prob_deriv[i]);
    EXPECT_FLOAT_EQ(b[i].adj(), dens_deriv[i]);
  }
};

TEST(AgradRevMatrix, log_mix_vals) {
  AVEC avec_prob(4), avec_dens(4);
  vector_v vecv_prob(4), vecv_dens(4);
  row_vector_v row_vecv_prob(4), row_vecv_dens(4);
  std::vector<var> stvec_prob(4), stvec_dens(4);

  val_rev_test(avec_prob, avec_dens);
  val_rev_test(avec_prob, vecv_dens);
  val_rev_test(avec_prob, row_vecv_dens);
  val_rev_test(avec_prob, stvec_dens);

  val_rev_test(vecv_prob, avec_dens);
  val_rev_test(vecv_prob, vecv_dens);
  val_rev_test(vecv_prob, row_vecv_dens);
  val_rev_test(vecv_prob, stvec_dens);

  val_rev_test(row_vecv_prob, avec_dens);
  val_rev_test(row_vecv_prob, vecv_dens);
  val_rev_test(row_vecv_prob, row_vecv_dens);
  val_rev_test(row_vecv_prob, stvec_dens);

  val_rev_test(stvec_prob, avec_dens);
  val_rev_test(stvec_prob, vecv_dens);
  val_rev_test(stvec_prob, row_vecv_dens);
  val_rev_test(stvec_prob, stvec_dens);
}
