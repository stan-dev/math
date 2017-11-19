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

template<typename T_prob, typename T_dens>
void mix_fv_fv(const T_prob& theta, const T_dens& lambda) {
  stan::math::vector_d prob_deriv(4);
  prob_deriv << 2.3610604993, 0.8685856170, 0.3195347914, 0.1175502804;

  stan::math::vector_d dens_deriv(4);
  dens_deriv << 0.3541590748, 0.6080099319, 0.0319534791, 0.0058775140;

  fvar<var> out = log_mix(theta, lambda);

  EXPECT_FLOAT_EQ(out.val_.val(), -1.85911088);
  EXPECT_FLOAT_EQ(out.d_.val(), 4.66673118);
  out.d_.grad();

  for (int i = 0; i < 4; ++i) {
    EXPECT_FLOAT_EQ(theta[i].d_.adj(), prob_deriv[i]);
    EXPECT_FLOAT_EQ(lambda[i].d_.adj(), dens_deriv[i]);
  }
}

template<typename T_prob, typename T_dens>
void mix_fv_d(const T_prob& theta, const T_dens& lambda) {
  stan::math::vector_d prob_deriv(4);
  prob_deriv << 2.3610604993, 0.8685856170, 0.3195347914, 0.1175502804;

  fvar<var> out = log_mix(theta, lambda);

  EXPECT_FLOAT_EQ(out.val_.val(), -1.85911088);
  EXPECT_FLOAT_EQ(out.d_.val(), 3.66673118);
  out.d_.grad();

  for (int i = 0; i < 4; ++i) {
    EXPECT_FLOAT_EQ(theta[i].d_.adj(), prob_deriv[i]);
  }
}

template<typename T_prob, typename T_dens>
void mix_d_fv(const T_prob& theta, const T_dens& lambda) {
  stan::math::vector_d dens_deriv(4);
  dens_deriv << 0.3541590748, 0.6080099319, 0.0319534791, 0.0058775140;

  fvar<var> out = log_mix(theta, lambda);

  EXPECT_FLOAT_EQ(out.val_.val(), -1.85911088);
  EXPECT_FLOAT_EQ(out.d_.val(), 1.0);
  out.d_.grad();

  for (int i = 0; i < 4; ++i) {
    EXPECT_FLOAT_EQ(lambda[i].d_.adj(), dens_deriv[i]);
  }
}

template<typename T_prob, typename T_dens>
void mix_ffv_ffv(const T_prob& theta, const T_dens& lambda) {
  stan::math::vector_d prob_deriv(4);
  prob_deriv << 2.3610604993, 0.8685856170, 0.3195347914, 0.1175502804;

  stan::math::vector_d dens_deriv(4);
  dens_deriv << 0.3541590748, 0.6080099319, 0.0319534791, 0.0058775140;

  fvar<fvar<var> > out = log_mix(theta, lambda);

  EXPECT_FLOAT_EQ(out.val_.val_.val(), -1.85911088);
  EXPECT_FLOAT_EQ(out.d_.val_.val(), 4.66673118);
  out.d_.val_.grad();

  for (int i = 0; i < 4; ++i) {
    EXPECT_FLOAT_EQ(theta[i].d_.val_.adj(), prob_deriv[i]);
    EXPECT_FLOAT_EQ(lambda[i].d_.val_.adj(), dens_deriv[i]);
  }
}

template<typename T_prob, typename T_dens>
void mix_ffv_d(const T_prob& theta, const T_dens& lambda) {
  stan::math::vector_d prob_deriv(4);
  prob_deriv << 2.3610604993, 0.8685856170, 0.3195347914, 0.1175502804;

  fvar<fvar<var> > out = log_mix(theta, lambda);

  EXPECT_FLOAT_EQ(out.val_.val_.val(), -1.85911088);
  EXPECT_FLOAT_EQ(out.d_.val_.val(), 3.66673118);
  out.d_.val_.grad();

  for (int i = 0; i < 4; ++i) {
    EXPECT_FLOAT_EQ(theta[i].d_.val_.adj(), prob_deriv[i]);
  }
}

template<typename T_prob, typename T_dens>
void mix_d_ffv(const T_prob& theta, const T_dens& lambda) {
  stan::math::vector_d dens_deriv(4);
  dens_deriv << 0.3541590748, 0.6080099319, 0.0319534791, 0.0058775140;

  fvar<fvar<var> > out = log_mix(theta, lambda);

  EXPECT_FLOAT_EQ(out.val_.val_.val(), -1.85911088);
  EXPECT_FLOAT_EQ(out.d_.val_.val(), 1.0);
  out.d_.val_.grad();

  for (int i = 0; i < 4; ++i) {
    EXPECT_FLOAT_EQ(lambda[i].d_.val_.adj(), dens_deriv[i]);
  }
}


TEST(AgradMixMatrixLogMix, fv_fv) {
  auto prob_fv = [](auto a) {
    a[0].val_ = 0.15;
    a[1].val_ = 0.70;
    a[2].val_ = 0.10;
    a[3].val_ = 0.05;
    a[0].d_ = 1.0;
    a[1].d_ = 1.0;
    a[2].d_ = 1.0;
    a[3].d_ = 1.0;

    return a;
  };
  auto dens_fv = [](auto b) {
    b[0].val_ = -1.0;
    b[1].val_ = -2.0;
    b[2].val_ = -3.0;
    b[3].val_ = -4.0;
    b[0].d_ = 1.0;
    b[1].d_ = 1.0;
    b[2].d_ = 1.0;
    b[3].d_ = 1.0;

    return b;
  };

  vector_fv vecfv_prob(4);
  vector_fv vecfv_dens(4);
  row_vector_fv row_vecfv_prob(4);
  row_vector_fv row_vecfv_dens(4);
  std::vector<fvar<var> > std_vecfv_prob(4);
  std::vector<fvar<var> > std_vecfv_dens(4);

  mix_fv_fv(prob_fv(vecfv_prob), dens_fv(vecfv_dens));
  mix_fv_fv(prob_fv(vecfv_prob), dens_fv(row_vecfv_dens));
  mix_fv_fv(prob_fv(vecfv_prob), dens_fv(std_vecfv_dens));

  mix_fv_fv(prob_fv(row_vecfv_prob), dens_fv(vecfv_dens));
  mix_fv_fv(prob_fv(row_vecfv_prob), dens_fv(row_vecfv_dens));
  mix_fv_fv(prob_fv(row_vecfv_prob), dens_fv(std_vecfv_dens));

  mix_fv_fv(prob_fv(std_vecfv_prob), dens_fv(vecfv_dens));
  mix_fv_fv(prob_fv(std_vecfv_prob), dens_fv(row_vecfv_dens));
  mix_fv_fv(prob_fv(std_vecfv_prob), dens_fv(std_vecfv_dens));
}

TEST(AgradMixMatrixLogMix, fv_d) {
  auto prob_fv = [](auto a) {
    a[0].val_ = 0.15;
    a[1].val_ = 0.70;
    a[2].val_ = 0.10;
    a[3].val_ = 0.05;
    a[0].d_ = 1.0;
    a[1].d_ = 1.0;
    a[2].d_ = 1.0;
    a[3].d_ = 1.0;

    return a;
  };

  auto dens_d = [](auto b) {
    b[0] = -1.0;
    b[1] = -2.0;
    b[2] = -3.0;
    b[3] = -4.0;

    return b;
  };

  vector_fv vecfv_prob(4);
  vector_d vecd_dens(4);
  row_vector_fv row_vecfv_prob(4);
  row_vector_d row_vecd_dens(4);
  std::vector<fvar<var> > std_vecfv_prob(4);
  std::vector<double> std_vecd_dens(4);

  mix_fv_d(prob_fv(vecfv_prob), dens_d(vecd_dens));
  mix_fv_d(prob_fv(vecfv_prob), dens_d(row_vecd_dens));
  mix_fv_d(prob_fv(vecfv_prob), dens_d(std_vecd_dens));

  mix_fv_d(prob_fv(row_vecfv_prob), dens_d(vecd_dens));
  mix_fv_d(prob_fv(row_vecfv_prob), dens_d(row_vecd_dens));
  mix_fv_d(prob_fv(row_vecfv_prob), dens_d(std_vecd_dens));

  mix_fv_d(prob_fv(std_vecfv_prob), dens_d(vecd_dens));
  mix_fv_d(prob_fv(std_vecfv_prob), dens_d(row_vecd_dens));
  mix_fv_d(prob_fv(std_vecfv_prob), dens_d(std_vecd_dens));
}

TEST(AgradMixMatrixLogMix, d_fv) {
  auto prob_d = [](auto a) {
    a[0] = 0.15;
    a[1] = 0.70;
    a[2] = 0.10;
    a[3] = 0.05;

    return a;
  };

  auto dens_fv = [](auto b) {
    b[0].val_ = -1.0;
    b[1].val_ = -2.0;
    b[2].val_ = -3.0;
    b[3].val_ = -4.0;
    b[0].d_ = 1.0;
    b[1].d_ = 1.0;
    b[2].d_ = 1.0;
    b[3].d_ = 1.0;

    return b;
  };

  vector_d vecd_prob(4);
  vector_fv vecfv_dens(4);
  row_vector_d row_vecd_prob(4);
  row_vector_fv row_vecfv_dens(4);
  std::vector<double> std_vecd_prob(4);
  std::vector<fvar<var> > std_vecfv_dens(4);

  mix_d_fv(prob_d(vecd_prob), dens_fv(vecfv_dens));
  mix_d_fv(prob_d(vecd_prob), dens_fv(row_vecfv_dens));
  mix_d_fv(prob_d(vecd_prob), dens_fv(std_vecfv_dens));

  mix_d_fv(prob_d(row_vecd_prob), dens_fv(vecfv_dens));
  mix_d_fv(prob_d(row_vecd_prob), dens_fv(row_vecfv_dens));
  mix_d_fv(prob_d(row_vecd_prob), dens_fv(std_vecfv_dens));

  mix_d_fv(prob_d(std_vecd_prob), dens_fv(vecfv_dens));
  mix_d_fv(prob_d(std_vecd_prob), dens_fv(row_vecfv_dens));
  mix_d_fv(prob_d(std_vecd_prob), dens_fv(std_vecfv_dens));
}

TEST(AgradMixMatrixLogMix, ffv_ffv) {
  auto prob_ffv = [](auto a) {
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

    return a;
  };
  auto dens_ffv = [](auto b) {
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

    return b;
  };

  vector_ffv vecffv_prob(4);
  vector_ffv vecffv_dens(4);
  row_vector_ffv row_vecffv_prob(4);
  row_vector_ffv row_vecffv_dens(4);
  std::vector<fvar<fvar<var> > > std_vecffv_prob(4);
  std::vector<fvar<fvar<var> > > std_vecffv_dens(4);

  mix_ffv_ffv(prob_ffv(vecffv_prob), dens_ffv(vecffv_dens));
  mix_ffv_ffv(prob_ffv(vecffv_prob), dens_ffv(row_vecffv_dens));
  mix_ffv_ffv(prob_ffv(vecffv_prob), dens_ffv(std_vecffv_dens));

  mix_ffv_ffv(prob_ffv(row_vecffv_prob), dens_ffv(vecffv_dens));
  mix_ffv_ffv(prob_ffv(row_vecffv_prob), dens_ffv(row_vecffv_dens));
  mix_ffv_ffv(prob_ffv(row_vecffv_prob), dens_ffv(std_vecffv_dens));

  mix_ffv_ffv(prob_ffv(std_vecffv_prob), dens_ffv(vecffv_dens));
  mix_ffv_ffv(prob_ffv(std_vecffv_prob), dens_ffv(row_vecffv_dens));
  mix_ffv_ffv(prob_ffv(std_vecffv_prob), dens_ffv(std_vecffv_dens));
}

TEST(AgradMixMatrixLogMix, ffv_d) {
  auto prob_ffv = [](auto a) {
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

    return a;
  };

  auto dens_d = [](auto b) {
    b[0] = -1.0;
    b[1] = -2.0;
    b[2] = -3.0;
    b[3] = -4.0;

    return b;
  };

  vector_ffv vecffv_prob(4);
  vector_d vecd_dens(4);
  row_vector_ffv row_vecffv_prob(4);
  row_vector_d row_vecd_dens(4);
  std::vector<fvar<fvar<var> > > std_vecffv_prob(4);
  std::vector<double> std_vecd_dens(4);

  mix_ffv_d(prob_ffv(vecffv_prob), dens_d(vecd_dens));
  mix_ffv_d(prob_ffv(vecffv_prob), dens_d(row_vecd_dens));
  mix_ffv_d(prob_ffv(vecffv_prob), dens_d(std_vecd_dens));

  mix_ffv_d(prob_ffv(row_vecffv_prob), dens_d(vecd_dens));
  mix_ffv_d(prob_ffv(row_vecffv_prob), dens_d(row_vecd_dens));
  mix_ffv_d(prob_ffv(row_vecffv_prob), dens_d(std_vecd_dens));

  mix_ffv_d(prob_ffv(std_vecffv_prob), dens_d(vecd_dens));
  mix_ffv_d(prob_ffv(std_vecffv_prob), dens_d(row_vecd_dens));
  mix_ffv_d(prob_ffv(std_vecffv_prob), dens_d(std_vecd_dens));
}

TEST(AgradMixMatrixLogMix, d_ffv) {
  auto prob_d = [](auto a) {
    a[0] = 0.15;
    a[1] = 0.70;
    a[2] = 0.10;
    a[3] = 0.05;

    return a;
  };

  auto dens_ffv = [](auto b) {
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

    return b;
  };

  vector_d vecd_prob(4);
  vector_ffv vecffv_dens(4);
  row_vector_d row_vecd_prob(4);
  row_vector_ffv row_vecffv_dens(4);
  std::vector<double> std_vecd_prob(4);
  std::vector<fvar<fvar<var> > > std_vecffv_dens(4);

  mix_d_ffv(prob_d(vecd_prob), dens_ffv(vecffv_dens));
  mix_d_ffv(prob_d(vecd_prob), dens_ffv(row_vecffv_dens));
  mix_d_ffv(prob_d(vecd_prob), dens_ffv(std_vecffv_dens));

  mix_d_ffv(prob_d(row_vecd_prob), dens_ffv(vecffv_dens));
  mix_d_ffv(prob_d(row_vecd_prob), dens_ffv(row_vecffv_dens));
  mix_d_ffv(prob_d(row_vecd_prob), dens_ffv(std_vecffv_dens));

  mix_d_ffv(prob_d(std_vecd_prob), dens_ffv(vecffv_dens));
  mix_d_ffv(prob_d(std_vecd_prob), dens_ffv(row_vecffv_dens));
  mix_d_ffv(prob_d(std_vecd_prob), dens_ffv(std_vecffv_dens));
}
