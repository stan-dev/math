#include <stan/math/mix/mat.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <gtest/gtest.h>

using stan::math::fvar;
using stan::math::var;
using stan::math::log_mix;

  template<typename T1>
  void check_prob_d(const T1& prob){
    stan::math::vector_d prob_d(4);
    prob_d << 2.3610604993, 0.8685856170, 0.3195347914, 0.1175502804;

    for (int i = 0; i < 4; ++i) {
      EXPECT_FLOAT_EQ(prob[i].d_.adj(), prob_d[i]);
    }
  }

  template<typename T2>
  void check_dens_d(const T2& dens){
    stan::math::vector_d dens_d(4);
    dens_d << 0.3541590748, 0.6080099319, 0.0319534791, 0.0058775140;

    for (int i = 0; i < 4; ++i) {
      EXPECT_FLOAT_EQ(dens[i].d_.adj(), dens_d[i]);
    }
  }
   

TEST(AgradMixMatrixLogMix, vec_fv_vec_fv) {
  using stan::math::vector_fv;

  vector_fv dens(4);
  dens << -1.0, -2.0, -3.0, -4.0;
  dens(0).d_ = 1.0;
  dens(1).d_ = 1.0;
  dens(2).d_ = 1.0;
  dens(3).d_ = 1.0;

  vector_fv prob(4);
  prob << 0.15, 0.70, 0.10, 0.05;
  prob(0).d_ = 1.0;
  prob(1).d_ = 1.0;
  prob(2).d_ = 1.0;
  prob(3).d_ = 1.0;

  fvar<var> out = log_mix(prob, dens);

  EXPECT_FLOAT_EQ(out.val_.val(), -1.85911088);
  EXPECT_FLOAT_EQ(out.d_.val(), 4.66673118);

  out.d_.grad();

  check_prob_d(prob);
  check_dens_d(dens);
}



TEST(AgradMixMatrixLogMix, rowvec_fv_vec_fv) {
  using stan::math::vector_fv;
  using stan::math::row_vector_fv;

  row_vector_fv dens(4);
  dens << -1.0, -2.0, -3.0, -4.0;
  dens(0).d_ = 1.0;
  dens(1).d_ = 1.0;
  dens(2).d_ = 1.0;
  dens(3).d_ = 1.0;

  vector_fv prob(4);
  prob << 0.15, 0.70, 0.10, 0.05;
  prob(0).d_ = 1.0;
  prob(1).d_ = 1.0;
  prob(2).d_ = 1.0;
  prob(3).d_ = 1.0;

  fvar<var> out = log_mix(prob, dens);

  EXPECT_FLOAT_EQ(out.val_.val(), -1.85911088);
  EXPECT_FLOAT_EQ(out.d_.val(), 4.66673118);

  out.d_.grad();

  check_prob_d(prob);
  check_dens_d(dens);
}

TEST(AgradMixMatrixLogMix, vec_fv_rowvec_fv) {
  using stan::math::vector_fv;
  using stan::math::row_vector_fv;

  vector_fv dens(4);
  dens << -1.0, -2.0, -3.0, -4.0;
  dens(0).d_ = 1.0;
  dens(1).d_ = 1.0;
  dens(2).d_ = 1.0;
  dens(3).d_ = 1.0;

  row_vector_fv prob(4);
  prob << 0.15, 0.70, 0.10, 0.05;
  prob(0).d_ = 1.0;
  prob(1).d_ = 1.0;
  prob(2).d_ = 1.0;
  prob(3).d_ = 1.0;

  fvar<var> out = log_mix(prob, dens);

  EXPECT_FLOAT_EQ(out.val_.val(), -1.85911088);
  EXPECT_FLOAT_EQ(out.d_.val(), 4.66673118);

  out.d_.grad();

  check_prob_d(prob);
  check_dens_d(dens);
}