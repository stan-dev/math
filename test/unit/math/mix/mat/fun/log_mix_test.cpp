#include <stan/math/mix/mat.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <gtest/gtest.h>
#include <vector>

using stan::math::fvar;
using stan::math::var;
using stan::math::log_mix;
using stan::math::vector_fv;
using stan::math::vector_ffv;
using stan::math::vector_d;
using stan::math::row_vector_fv;
using stan::math::row_vector_ffv;
using stan::math::row_vector_d;

TEST(AgradMixMatrixLogMix, fv_fv) {
  auto mix_fv_fv = [](auto a, auto b) {
    a[0].val_ = 0.15;
    a[1].val_ = 0.70;
    a[2].val_ = 0.10;
    a[3].val_ = 0.05;
    a[0].d_ = 1.0;
    a[1].d_ = 1.0;
    a[2].d_ = 1.0;
    a[3].d_ = 1.0;

    b[0].val_ = -1.0;
    b[1].val_ = -2.0;
    b[2].val_ = -3.0;
    b[3].val_ = -4.0;
    b[0].d_ = 1.0;
    b[1].d_ = 1.0;
    b[2].d_ = 1.0;
    b[3].d_ = 1.0;

    fvar<var> out = log_mix(a, b);

    EXPECT_FLOAT_EQ(out.val_.val(), -1.85911088);
    EXPECT_FLOAT_EQ(out.d_.val(), 4.66673118);
    out.d_.grad();

    vector_d prob_deriv(4);
    prob_deriv << 2.3610604993, 0.8685856170, 0.3195347914, 0.1175502804;
    vector_d dens_deriv(4);
    dens_deriv << 0.3541590748, 0.6080099319, 0.0319534791, 0.0058775140;

    for (int i = 0; i < 4; ++i) {
      EXPECT_FLOAT_EQ(a[i].d_.adj(), prob_deriv[i]);
      EXPECT_FLOAT_EQ(b[i].d_.adj(), dens_deriv[i]);
    }
  };

  vector_fv vecfv_prob(4);
  vector_fv vecfv_dens(4);
  row_vector_fv row_vecfv_prob(4);
  row_vector_fv row_vecfv_dens(4);
  std::vector<fvar<var> > std_vecfv_prob(4);
  std::vector<fvar<var> > std_vecfv_dens(4);

  mix_fv_fv(vecfv_prob, vecfv_dens);
  mix_fv_fv(vecfv_prob, row_vecfv_dens);
  mix_fv_fv(vecfv_prob, std_vecfv_dens);

  mix_fv_fv(row_vecfv_prob, vecfv_dens);
  mix_fv_fv(row_vecfv_prob, row_vecfv_dens);
  mix_fv_fv(row_vecfv_prob, std_vecfv_dens);

  mix_fv_fv(std_vecfv_prob, vecfv_dens);
  mix_fv_fv(std_vecfv_prob, row_vecfv_dens);
  mix_fv_fv(std_vecfv_prob, std_vecfv_dens);
}

TEST(AgradMixMatrixLogMix, fv_d) {
  auto mix_fv_d = [](auto a, auto b) {
    a[0].val_ = 0.15;
    a[1].val_ = 0.70;
    a[2].val_ = 0.10;
    a[3].val_ = 0.05;
    a[0].d_ = 1.0;
    a[1].d_ = 1.0;
    a[2].d_ = 1.0;
    a[3].d_ = 1.0;

    b[0] = -1.0;
    b[1] = -2.0;
    b[2] = -3.0;
    b[3] = -4.0;

    fvar<var> out = log_mix(a, b);

    EXPECT_FLOAT_EQ(out.val_.val(), -1.85911088);
    EXPECT_FLOAT_EQ(out.d_.val(), 3.66673118);
    out.d_.grad();

    vector_d prob_deriv(4);
    prob_deriv << 2.3610604993, 0.8685856170, 0.3195347914, 0.1175502804;

    for (int i = 0; i < 4; ++i) {
      EXPECT_FLOAT_EQ(a[i].d_.adj(), prob_deriv[i]);
    }
  };

  vector_fv vecfv_prob(4);
  vector_d vecd_dens(4);
  row_vector_fv row_vecfv_prob(4);
  row_vector_d row_vecd_dens(4);
  std::vector<fvar<var> > std_vecfv_prob(4);
  std::vector<double> std_vecd_dens(4);

  mix_fv_d(vecfv_prob, vecd_dens);
  mix_fv_d(vecfv_prob, row_vecd_dens);
  mix_fv_d(vecfv_prob, std_vecd_dens);

  mix_fv_d(row_vecfv_prob, vecd_dens);
  mix_fv_d(row_vecfv_prob, row_vecd_dens);
  mix_fv_d(row_vecfv_prob, std_vecd_dens);

  mix_fv_d(std_vecfv_prob, vecd_dens);
  mix_fv_d(std_vecfv_prob, row_vecd_dens);
  mix_fv_d(std_vecfv_prob, std_vecd_dens);
}


TEST(AgradMixMatrixLogMix, d_fv) {
  auto mix_d_fv = [](auto a, auto b) {
    a[0] = 0.15;
    a[1] = 0.70;
    a[2] = 0.10;
    a[3] = 0.05;

    b[0].val_ = -1.0;
    b[1].val_ = -2.0;
    b[2].val_ = -3.0;
    b[3].val_ = -4.0;
    b[0].d_ = 1.0;
    b[1].d_ = 1.0;
    b[2].d_ = 1.0;
    b[3].d_ = 1.0;

    fvar<var> out = log_mix(a, b);

    EXPECT_FLOAT_EQ(out.val_.val(), -1.85911088);
    EXPECT_FLOAT_EQ(out.d_.val(), 1.0);
    out.d_.grad();

    vector_d dens_deriv(4);
    dens_deriv << 0.3541590748, 0.6080099319, 0.0319534791, 0.0058775140;

    for (int i = 0; i < 4; ++i) {
      EXPECT_FLOAT_EQ(b[i].d_.adj(), dens_deriv[i]);
    }
  };


  vector_d vecd_prob(4);
  vector_fv vecfv_dens(4);
  row_vector_d row_vecd_prob(4);
  row_vector_fv row_vecfv_dens(4);
  std::vector<double> std_vecd_prob(4);
  std::vector<fvar<var> > std_vecfv_dens(4);

  mix_d_fv(vecd_prob, vecfv_dens);
  mix_d_fv(vecd_prob, row_vecfv_dens);
  mix_d_fv(vecd_prob, std_vecfv_dens);

  mix_d_fv(row_vecd_prob, vecfv_dens);
  mix_d_fv(row_vecd_prob, row_vecfv_dens);
  mix_d_fv(row_vecd_prob, std_vecfv_dens);

  mix_d_fv(std_vecd_prob, vecfv_dens);
  mix_d_fv(std_vecd_prob, row_vecfv_dens);
  mix_d_fv(std_vecd_prob, std_vecfv_dens);
}

TEST(AgradMixMatrixLogMix, ffv_ffv) {
  auto mix_ffv_ffv = [](auto a, auto b) {
    a[0].val_ = 0.15;
    a[1].val_ = 0.70;
    a[2].val_ = 0.10;
    a[3].val_ = 0.05;
    a[0].d_ = 1.0;
    a[1].d_ = 1.0;
    a[2].d_ = 1.0;
    a[3].d_ = 1.0;
    a[0].val_.d_ = 1.0;
    a[1].val_.d_ = 1.0;
    a[2].val_.d_ = 1.0;
    a[3].val_.d_ = 1.0;

    b[0].val_ = -1.0;
    b[1].val_ = -2.0;
    b[2].val_ = -3.0;
    b[3].val_ = -4.0;
    b[0].d_ = 1.0;
    b[1].d_ = 1.0;
    b[2].d_ = 1.0;
    b[3].d_ = 1.0;
    b[0].val_.d_ = 1.0;
    b[1].val_.d_ = 1.0;
    b[2].val_.d_ = 1.0;
    b[3].val_.d_ = 1.0;

    fvar<fvar<var> > out = log_mix(a, b);

    EXPECT_FLOAT_EQ(out.val_.val_.val(), -1.85911088);
    EXPECT_FLOAT_EQ(out.d_.val_.val(), 4.66673118);
    out.d_.val_.grad();

    stan::math::vector_d prob_deriv(4);
    prob_deriv << 2.3610604993, 0.8685856170, 0.3195347914, 0.1175502804;
    stan::math::vector_d dens_deriv(4);
    dens_deriv << 0.3541590748, 0.6080099319, 0.0319534791, 0.0058775140;

    for (int i = 0; i < 4; ++i) {
      EXPECT_FLOAT_EQ(a[i].d_.val_.adj(), prob_deriv[i]);
      EXPECT_FLOAT_EQ(b[i].d_.val_.adj(), dens_deriv[i]);
    }
  };

  vector_ffv vecffv_prob(4);
  vector_ffv vecffv_dens(4);
  row_vector_ffv row_vecffv_prob(4);
  row_vector_ffv row_vecffv_dens(4);
  std::vector<fvar<fvar<var> > > std_vecffv_prob(4);
  std::vector<fvar<fvar<var> > > std_vecffv_dens(4);

  mix_ffv_ffv(vecffv_prob, vecffv_dens);
  mix_ffv_ffv(vecffv_prob, row_vecffv_dens);
  mix_ffv_ffv(vecffv_prob, std_vecffv_dens);

  mix_ffv_ffv(row_vecffv_prob, vecffv_dens);
  mix_ffv_ffv(row_vecffv_prob, row_vecffv_dens);
  mix_ffv_ffv(row_vecffv_prob, std_vecffv_dens);

  mix_ffv_ffv(std_vecffv_prob, vecffv_dens);
  mix_ffv_ffv(std_vecffv_prob, row_vecffv_dens);
  mix_ffv_ffv(std_vecffv_prob, std_vecffv_dens);
}

TEST(AgradMixMatrixLogMix, ffv_d) {
  auto mix_ffv_d = [](auto a, auto b) {
    a[0].val_ = 0.15;
    a[1].val_ = 0.70;
    a[2].val_ = 0.10;
    a[3].val_ = 0.05;
    a[0].d_ = 1.0;
    a[1].d_ = 1.0;
    a[2].d_ = 1.0;
    a[3].d_ = 1.0;
    a[0].val_.d_ = 1.0;
    a[1].val_.d_ = 1.0;
    a[2].val_.d_ = 1.0;
    a[3].val_.d_ = 1.0;

    b[0] = -1.0;
    b[1] = -2.0;
    b[2] = -3.0;
    b[3] = -4.0;

    fvar<fvar<var> > out = log_mix(a, b);

    EXPECT_FLOAT_EQ(out.val_.val_.val(), -1.85911088);
    EXPECT_FLOAT_EQ(out.d_.val_.val(), 3.66673118);
    out.d_.val_.grad();

    stan::math::vector_d prob_deriv(4);
    prob_deriv << 2.3610604993, 0.8685856170, 0.3195347914, 0.1175502804;

    for (int i = 0; i < 4; ++i) {
      EXPECT_FLOAT_EQ(a[i].d_.val_.adj(), prob_deriv[i]);
    }
  };

  vector_ffv vecffv_prob(4);
  vector_d vecd_dens(4);
  row_vector_ffv row_vecffv_prob(4);
  row_vector_d row_vecd_dens(4);
  std::vector<fvar<fvar<var> > > std_vecffv_prob(4);
  std::vector<double> std_vecd_dens(4);

  mix_ffv_d(vecffv_prob, vecd_dens);
  mix_ffv_d(vecffv_prob, row_vecd_dens);
  mix_ffv_d(vecffv_prob, std_vecd_dens);

  mix_ffv_d(row_vecffv_prob, vecd_dens);
  mix_ffv_d(row_vecffv_prob, row_vecd_dens);
  mix_ffv_d(row_vecffv_prob, std_vecd_dens);

  mix_ffv_d(std_vecffv_prob, vecd_dens);
  mix_ffv_d(std_vecffv_prob, row_vecd_dens);
  mix_ffv_d(std_vecffv_prob, std_vecd_dens);
}

TEST(AgradMixMatrixLogMix, d_ffv) {
  auto mix_d_ffv = [](auto a, auto b) {
    a[0] = 0.15;
    a[1] = 0.70;
    a[2] = 0.10;
    a[3] = 0.05;

    b[0].val_ = -1.0;
    b[1].val_ = -2.0;
    b[2].val_ = -3.0;
    b[3].val_ = -4.0;
    b[0].d_ = 1.0;
    b[1].d_ = 1.0;
    b[2].d_ = 1.0;
    b[3].d_ = 1.0;
    b[0].val_.d_ = 1.0;
    b[1].val_.d_ = 1.0;
    b[2].val_.d_ = 1.0;
    b[3].val_.d_ = 1.0;

    fvar<fvar<var> > out = log_mix(a, b);

    EXPECT_FLOAT_EQ(out.val_.val_.val(), -1.85911088);
    EXPECT_FLOAT_EQ(out.d_.val_.val(), 1.0);
    out.d_.val_.grad();

    stan::math::vector_d dens_deriv(4);
    dens_deriv << 0.3541590748, 0.6080099319, 0.0319534791, 0.0058775140;

    for (int i = 0; i < 4; ++i) {
      EXPECT_FLOAT_EQ(b[i].d_.val_.adj(), dens_deriv[i]);
    }
  };

  vector_d vecd_prob(4);
  vector_ffv vecffv_dens(4);
  row_vector_d row_vecd_prob(4);
  row_vector_ffv row_vecffv_dens(4);
  std::vector<double> std_vecd_prob(4);
  std::vector<fvar<fvar<var> > > std_vecffv_dens(4);

  mix_d_ffv(vecd_prob, vecffv_dens);
  mix_d_ffv(vecd_prob, row_vecffv_dens);
  mix_d_ffv(vecd_prob, std_vecffv_dens);

  mix_d_ffv(row_vecd_prob, vecffv_dens);
  mix_d_ffv(row_vecd_prob, row_vecffv_dens);
  mix_d_ffv(row_vecd_prob, std_vecffv_dens);

  mix_d_ffv(std_vecd_prob, vecffv_dens);
  mix_d_ffv(std_vecd_prob, row_vecffv_dens);
  mix_d_ffv(std_vecd_prob, std_vecffv_dens);
}
