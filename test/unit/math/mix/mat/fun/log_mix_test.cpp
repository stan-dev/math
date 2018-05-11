#include <stan/math/mix/mat.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <gtest/gtest.h>
#include <vector>

using stan::math::fvar;
using stan::math::log_mix;
using stan::math::row_vector_d;
using stan::math::row_vector_ffv;
using stan::math::row_vector_fv;
using stan::math::var;
using stan::math::vector_d;
using stan::math::vector_ffv;
using stan::math::vector_fv;

template <typename T_a, typename T_b>
void fv_fv_test(T_a a, T_b b) {
  a[0].val_ = 0.514;
  a[1].val_ = 0.284;
  a[2].val_ = 0.112;
  a[3].val_ = 0.090;
  a[0].d_ = 1.0;
  a[1].d_ = 1.0;
  a[2].d_ = 1.0;
  a[3].d_ = 1.0;

  b[0].val_ = -3.581;
  b[1].val_ = -8.114;
  b[2].val_ = -11.215;
  b[3].val_ = -5.658;
  b[0].d_ = 1.0;
  b[1].d_ = 1.0;
  b[2].d_ = 1.0;
  b[3].d_ = 1.0;

  fvar<var> out = log_mix(a, b);

  EXPECT_FLOAT_EQ(out.val_.val(), -4.218931574);
  EXPECT_FLOAT_EQ(out.d_.val(), 3.150968236);
  out.d_.grad();

  vector_d prob_deriv(4);
  prob_deriv << 1.892562202198, 0.020341982471, 0.000915474153, 0.237148577149;
  vector_d dens_deriv(4);
  dens_deriv << 0.97277697193, 0.00577712302, 0.00010253311, 0.02134337194;

  for (int i = 0; i < 4; ++i) {
    EXPECT_FLOAT_EQ(a[i].d_.adj(), prob_deriv[i]);
    EXPECT_FLOAT_EQ(b[i].d_.adj(), dens_deriv[i]);
  }
}

template <typename T_a, typename T_b>
void fv_fv_vec_test(T_a a, T_b b1) {
  a[0].val_ = 0.514;
  a[1].val_ = 0.284;
  a[2].val_ = 0.112;
  a[3].val_ = 0.090;
  a[0].d_ = 1.0;
  a[1].d_ = 1.0;
  a[2].d_ = 1.0;
  a[3].d_ = 1.0;

  b1[0].val_ = -3.581;
  b1[1].val_ = -8.114;
  b1[2].val_ = -11.215;
  b1[3].val_ = -5.658;
  b1[0].d_ = 1.0;
  b1[1].d_ = 1.0;
  b1[2].d_ = 1.0;
  b1[3].d_ = 1.0;

  T_b b2(4);

  b2[0].val_ = -8.594;
  b2[1].val_ = -3.251;
  b2[2].val_ = -7.281;
  b2[3].val_ = -3.556;
  b2[0].d_ = 1.0;
  b2[1].d_ = 1.0;
  b2[2].d_ = 1.0;
  b2[3].d_ = 1.0;

  T_b b3(4);

  b3[0].val_ = -11.554;
  b3[1].val_ = -6.628;
  b3[2].val_ = -15.229;
  b3[3].val_ = -9.561;
  b3[0].d_ = 1.0;
  b3[1].d_ = 1.0;
  b3[2].d_ = 1.0;
  b3[3].d_ = 1.0;

  std::vector<T_b> c{b1, b2, b3};

  fvar<var> out = log_mix(a, c);

  EXPECT_FLOAT_EQ(out.val_.val(), -16.36331174);
  EXPECT_FLOAT_EQ(out.d_.val(), 13.73648479);
  out.d_.grad();

  vector_d prob_deriv(4);
  prob_deriv << 1.930840738468, 6.257235890773, 0.051642416014, 2.496765742829;
  vector_d dens1_deriv(4);
  dens1_deriv << 0.97277697193, 0.00577712302, 0.00010253311, 0.02134337194;
  vector_d dens2_deriv(4);
  dens2_deriv << 0.00692718708, 0.80047459791, 0.00561100267, 0.18698721234;
  vector_d dens3_deriv(4);
  dens3_deriv << 0.01274798056, 0.97080327205, 0.00007041482, 0.01637833257;

  for (int i = 0; i < 4; ++i) {
    EXPECT_FLOAT_EQ(a[i].d_.adj(), prob_deriv[i]);
    EXPECT_FLOAT_EQ(c[0][i].d_.adj(), dens1_deriv[i]);
    EXPECT_FLOAT_EQ(c[1][i].d_.adj(), dens2_deriv[i]);
    EXPECT_FLOAT_EQ(c[2][i].d_.adj(), dens3_deriv[i]);
  }
}

TEST(AgradMixMatrixLogMix, fv_fv) {
  vector_fv vecfv_prob(4);
  vector_fv vecfv_dens(4);
  row_vector_fv row_vecfv_prob(4);
  row_vector_fv row_vecfv_dens(4);
  std::vector<fvar<var> > std_vecfv_prob(4);
  std::vector<fvar<var> > std_vecfv_dens(4);

  fv_fv_test(vecfv_prob, vecfv_dens);
  fv_fv_test(vecfv_prob, row_vecfv_dens);
  fv_fv_test(vecfv_prob, std_vecfv_dens);

  fv_fv_vec_test(vecfv_prob, vecfv_dens);
  fv_fv_vec_test(vecfv_prob, row_vecfv_dens);
  fv_fv_vec_test(vecfv_prob, std_vecfv_dens);

  fv_fv_test(row_vecfv_prob, vecfv_dens);
  fv_fv_test(row_vecfv_prob, row_vecfv_dens);
  fv_fv_test(row_vecfv_prob, std_vecfv_dens);

  fv_fv_vec_test(row_vecfv_prob, vecfv_dens);
  fv_fv_vec_test(row_vecfv_prob, row_vecfv_dens);
  fv_fv_vec_test(row_vecfv_prob, std_vecfv_dens);

  fv_fv_test(std_vecfv_prob, vecfv_dens);
  fv_fv_test(std_vecfv_prob, row_vecfv_dens);
  fv_fv_test(std_vecfv_prob, std_vecfv_dens);

  fv_fv_vec_test(std_vecfv_prob, vecfv_dens);
  fv_fv_vec_test(std_vecfv_prob, row_vecfv_dens);
  fv_fv_vec_test(std_vecfv_prob, std_vecfv_dens);
}

template <typename T_a, typename T_b>
void fv_d_test(T_a a, T_b b) {
  a[0].val_ = 0.514;
  a[1].val_ = 0.284;
  a[2].val_ = 0.112;
  a[3].val_ = 0.090;
  a[0].d_ = 1.0;
  a[1].d_ = 1.0;
  a[2].d_ = 1.0;
  a[3].d_ = 1.0;

  b[0] = -3.581;
  b[1] = -8.114;
  b[2] = -11.215;
  b[3] = -5.658;

  fvar<var> out = log_mix(a, b);

  EXPECT_FLOAT_EQ(out.val_.val(), -4.2189315474);
  EXPECT_FLOAT_EQ(out.d_.val(), 2.150968235971);
  out.d_.grad();

  vector_d prob_deriv(4);
  prob_deriv << 1.892562202198, 0.020341982471, 0.000915474153, 0.237148577149;

  for (int i = 0; i < 4; ++i) {
    EXPECT_FLOAT_EQ(a[i].d_.adj(), prob_deriv[i]);
  }
}

template <typename T_a, typename T_b>
void fv_d_vec_test(T_a a, T_b b1) {
  a[0].val_ = 0.514;
  a[1].val_ = 0.284;
  a[2].val_ = 0.112;
  a[3].val_ = 0.090;
  a[0].d_ = 1.0;
  a[1].d_ = 1.0;
  a[2].d_ = 1.0;
  a[3].d_ = 1.0;

  b1[0] = -3.581;
  b1[1] = -8.114;
  b1[2] = -11.215;
  b1[3] = -5.658;

  T_b b2(4);

  b2[0] = -8.594;
  b2[1] = -3.251;
  b2[2] = -7.281;
  b2[3] = -3.556;

  T_b b3(4);

  b3[0] = -11.554;
  b3[1] = -6.628;
  b3[2] = -15.229;
  b3[3] = -9.561;

  std::vector<T_b> c{b1, b2, b3};

  fvar<var> out = log_mix(a, c);

  EXPECT_FLOAT_EQ(out.val_.val(), -16.36331174);
  EXPECT_FLOAT_EQ(out.d_.val(), 10.73648479);
  out.d_.grad();

  vector_d prob_deriv(4);
  prob_deriv << 1.930840738468, 6.257235890773, 0.051642416014, 2.496765742829;

  for (int i = 0; i < 4; ++i) {
    EXPECT_FLOAT_EQ(a[i].d_.adj(), prob_deriv[i]);
  }
}

TEST(AgradMixMatrixLogMix, fv_d) {
  vector_fv vecfv_prob(4);
  vector_d vecd_dens(4);
  row_vector_fv row_vecfv_prob(4);
  row_vector_d row_vecd_dens(4);
  std::vector<fvar<var> > std_vecfv_prob(4);
  std::vector<double> std_vecd_dens(4);

  fv_d_test(vecfv_prob, vecd_dens);
  fv_d_test(vecfv_prob, row_vecd_dens);
  fv_d_test(vecfv_prob, std_vecd_dens);

  fv_d_vec_test(vecfv_prob, vecd_dens);
  fv_d_vec_test(vecfv_prob, row_vecd_dens);
  fv_d_vec_test(vecfv_prob, std_vecd_dens);

  fv_d_test(row_vecfv_prob, vecd_dens);
  fv_d_test(row_vecfv_prob, row_vecd_dens);
  fv_d_test(row_vecfv_prob, std_vecd_dens);

  fv_d_vec_test(row_vecfv_prob, vecd_dens);
  fv_d_vec_test(row_vecfv_prob, row_vecd_dens);
  fv_d_vec_test(row_vecfv_prob, std_vecd_dens);

  fv_d_test(std_vecfv_prob, vecd_dens);
  fv_d_test(std_vecfv_prob, row_vecd_dens);
  fv_d_test(std_vecfv_prob, std_vecd_dens);

  fv_d_vec_test(std_vecfv_prob, vecd_dens);
  fv_d_vec_test(std_vecfv_prob, row_vecd_dens);
  fv_d_vec_test(std_vecfv_prob, std_vecd_dens);
}

template <typename T_a, typename T_b>
void d_fv_test(T_a a, T_b b) {
  a[0] = 0.514;
  a[1] = 0.284;
  a[2] = 0.112;
  a[3] = 0.090;

  b[0].val_ = -3.581;
  b[1].val_ = -8.114;
  b[2].val_ = -11.215;
  b[3].val_ = -5.658;
  b[0].d_ = 1.0;
  b[1].d_ = 1.0;
  b[2].d_ = 1.0;
  b[3].d_ = 1.0;

  fvar<var> out = log_mix(a, b);

  EXPECT_FLOAT_EQ(out.val_.val(), -4.218931574);
  EXPECT_FLOAT_EQ(out.d_.val(), 1.0);
  out.d_.grad();

  vector_d dens_deriv(4);
  dens_deriv << 0.97277697193, 0.00577712302, 0.00010253311, 0.02134337194;

  for (int i = 0; i < 4; ++i) {
    EXPECT_FLOAT_EQ(b[i].d_.adj(), dens_deriv[i]);
  }
}

template <typename T_a, typename T_b>
void d_fv_vec_test(T_a a, T_b b1) {
  a[0] = 0.514;
  a[1] = 0.284;
  a[2] = 0.112;
  a[3] = 0.090;

  b1[0].val_ = -3.581;
  b1[1].val_ = -8.114;
  b1[2].val_ = -11.215;
  b1[3].val_ = -5.658;
  b1[0].d_ = 1.0;
  b1[1].d_ = 1.0;
  b1[2].d_ = 1.0;
  b1[3].d_ = 1.0;

  T_b b2(4);

  b2[0].val_ = -8.594;
  b2[1].val_ = -3.251;
  b2[2].val_ = -7.281;
  b2[3].val_ = -3.556;
  b2[0].d_ = 1.0;
  b2[1].d_ = 1.0;
  b2[2].d_ = 1.0;
  b2[3].d_ = 1.0;

  T_b b3(4);

  b3[0].val_ = -11.554;
  b3[1].val_ = -6.628;
  b3[2].val_ = -15.229;
  b3[3].val_ = -9.561;
  b3[0].d_ = 1.0;
  b3[1].d_ = 1.0;
  b3[2].d_ = 1.0;
  b3[3].d_ = 1.0;

  std::vector<T_b> c{b1, b2, b3};

  fvar<var> out = log_mix(a, c);

  EXPECT_FLOAT_EQ(out.val_.val(), -16.36331174);
  EXPECT_FLOAT_EQ(out.d_.val(), 3.0);
  out.d_.grad();

  vector_d dens1_deriv(4);
  dens1_deriv << 0.97277697193, 0.00577712302, 0.00010253311, 0.02134337194;
  vector_d dens2_deriv(4);
  dens2_deriv << 0.00692718708, 0.80047459791, 0.00561100267, 0.18698721234;
  vector_d dens3_deriv(4);
  dens3_deriv << 0.01274798056, 0.97080327205, 0.00007041482, 0.01637833257;

  for (int i = 0; i < 4; ++i) {
    EXPECT_FLOAT_EQ(c[0][i].d_.adj(), dens1_deriv[i]);
    EXPECT_FLOAT_EQ(c[1][i].d_.adj(), dens2_deriv[i]);
    EXPECT_FLOAT_EQ(c[2][i].d_.adj(), dens3_deriv[i]);
  }
}

TEST(AgradMixMatrixLogMix, d_fv) {
  vector_d vecd_prob(4);
  vector_fv vecfv_dens(4);
  row_vector_d row_vecd_prob(4);
  row_vector_fv row_vecfv_dens(4);
  std::vector<double> std_vecd_prob(4);
  std::vector<fvar<var> > std_vecfv_dens(4);

  d_fv_test(vecd_prob, vecfv_dens);
  d_fv_test(vecd_prob, row_vecfv_dens);
  d_fv_test(vecd_prob, std_vecfv_dens);

  d_fv_vec_test(vecd_prob, vecfv_dens);
  d_fv_vec_test(vecd_prob, row_vecfv_dens);
  d_fv_vec_test(vecd_prob, std_vecfv_dens);

  d_fv_test(row_vecd_prob, vecfv_dens);
  d_fv_test(row_vecd_prob, row_vecfv_dens);
  d_fv_test(row_vecd_prob, std_vecfv_dens);

  d_fv_vec_test(row_vecd_prob, vecfv_dens);
  d_fv_vec_test(row_vecd_prob, row_vecfv_dens);
  d_fv_vec_test(row_vecd_prob, std_vecfv_dens);

  d_fv_test(std_vecd_prob, vecfv_dens);
  d_fv_test(std_vecd_prob, row_vecfv_dens);
  d_fv_test(std_vecd_prob, std_vecfv_dens);

  d_fv_vec_test(std_vecd_prob, vecfv_dens);
  d_fv_vec_test(std_vecd_prob, row_vecfv_dens);
  d_fv_vec_test(std_vecd_prob, std_vecfv_dens);
}

template <typename T_a, typename T_b>
void ffv_ffv_test(T_a a, T_b b) {
  a[0].val_ = 0.514;
  a[1].val_ = 0.284;
  a[2].val_ = 0.112;
  a[3].val_ = 0.090;
  a[0].d_ = 1.0;
  a[1].d_ = 1.0;
  a[2].d_ = 1.0;
  a[3].d_ = 1.0;
  a[0].val_.d_ = 1.0;
  a[1].val_.d_ = 1.0;
  a[2].val_.d_ = 1.0;
  a[3].val_.d_ = 1.0;

  b[0].val_ = -3.581;
  b[1].val_ = -8.114;
  b[2].val_ = -11.215;
  b[3].val_ = -5.658;
  b[0].d_ = 1.0;
  b[1].d_ = 1.0;
  b[2].d_ = 1.0;
  b[3].d_ = 1.0;
  b[0].val_.d_ = 1.0;
  b[1].val_.d_ = 1.0;
  b[2].val_.d_ = 1.0;
  b[3].val_.d_ = 1.0;

  fvar<fvar<var> > out = log_mix(a, b);

  EXPECT_FLOAT_EQ(out.val_.val_.val(), -4.218931574);
  EXPECT_FLOAT_EQ(out.d_.val_.val(), 3.150968236);
  out.d_.val_.grad();

  vector_d prob_deriv(4);
  prob_deriv << 1.892562202198, 0.020341982471, 0.000915474153, 0.237148577149;
  vector_d dens_deriv(4);
  dens_deriv << 0.97277697193, 0.00577712302, 0.00010253311, 0.02134337194;

  for (int i = 0; i < 4; ++i) {
    EXPECT_FLOAT_EQ(a[i].d_.val_.adj(), prob_deriv[i]);
    EXPECT_FLOAT_EQ(b[i].d_.val_.adj(), dens_deriv[i]);
  }
}

template <typename T_a, typename T_b>
void ffv_ffv_vec_test(T_a a, T_b b1) {
  a[0].val_ = 0.514;
  a[1].val_ = 0.284;
  a[2].val_ = 0.112;
  a[3].val_ = 0.090;
  a[0].d_ = 1.0;
  a[1].d_ = 1.0;
  a[2].d_ = 1.0;
  a[3].d_ = 1.0;
  a[0].val_.d_ = 1.0;
  a[1].val_.d_ = 1.0;
  a[2].val_.d_ = 1.0;
  a[3].val_.d_ = 1.0;

  b1[0].val_ = -3.581;
  b1[1].val_ = -8.114;
  b1[2].val_ = -11.215;
  b1[3].val_ = -5.658;

  b1[0].d_ = 1.0;
  b1[1].d_ = 1.0;
  b1[2].d_ = 1.0;
  b1[3].d_ = 1.0;
  b1[0].val_.d_ = 1.0;
  b1[1].val_.d_ = 1.0;
  b1[2].val_.d_ = 1.0;
  b1[3].val_.d_ = 1.0;

  T_b b2(4), b3(4);

  b2[0].val_ = -8.594;
  b2[1].val_ = -3.251;
  b2[2].val_ = -7.281;
  b2[3].val_ = -3.556;

  b2[0].d_ = 1.0;
  b2[1].d_ = 1.0;
  b2[2].d_ = 1.0;
  b2[3].d_ = 1.0;
  b2[0].val_.d_ = 1.0;
  b2[1].val_.d_ = 1.0;
  b2[2].val_.d_ = 1.0;
  b2[3].val_.d_ = 1.0;

  b3[0].val_ = -11.554;
  b3[1].val_ = -6.628;
  b3[2].val_ = -15.229;
  b3[3].val_ = -9.561;

  b3[0].d_ = 1.0;
  b3[1].d_ = 1.0;
  b3[2].d_ = 1.0;
  b3[3].d_ = 1.0;
  b3[0].val_.d_ = 1.0;
  b3[1].val_.d_ = 1.0;
  b3[2].val_.d_ = 1.0;
  b3[3].val_.d_ = 1.0;

  std::vector<T_b> c{b1, b2, b3};

  fvar<fvar<var> > out = log_mix(a, c);

  EXPECT_FLOAT_EQ(out.val_.val_.val(), -16.36331174);
  EXPECT_FLOAT_EQ(out.d_.val_.val(), 13.73648479);
  out.d_.val_.grad();

  vector_d prob_deriv(4);
  prob_deriv << 1.930840738468, 6.257235890773, 0.051642416014, 2.496765742829;
  vector_d dens1_deriv(4);
  dens1_deriv << 0.97277697193, 0.00577712302, 0.00010253311, 0.02134337194;
  vector_d dens2_deriv(4);
  dens2_deriv << 0.00692718708, 0.80047459791, 0.00561100267, 0.18698721234;
  vector_d dens3_deriv(4);
  dens3_deriv << 0.01274798056, 0.97080327205, 0.00007041482, 0.01637833257;

  for (int i = 0; i < 4; ++i) {
    EXPECT_FLOAT_EQ(a[i].d_.val_.adj(), prob_deriv[i]);
    EXPECT_FLOAT_EQ(c[0][i].d_.val_.adj(), dens1_deriv[i]);
    EXPECT_FLOAT_EQ(c[1][i].d_.val_.adj(), dens2_deriv[i]);
    EXPECT_FLOAT_EQ(c[2][i].d_.val_.adj(), dens3_deriv[i]);
  }
}

TEST(AgradMixMatrixLogMix, ffv_ffv) {
  vector_ffv vecffv_prob(4);
  vector_ffv vecffv_dens(4);
  row_vector_ffv row_vecffv_prob(4);
  row_vector_ffv row_vecffv_dens(4);
  std::vector<fvar<fvar<var> > > std_vecffv_prob(4);
  std::vector<fvar<fvar<var> > > std_vecffv_dens(4);

  ffv_ffv_test(vecffv_prob, vecffv_dens);
  ffv_ffv_test(vecffv_prob, row_vecffv_dens);
  ffv_ffv_test(vecffv_prob, std_vecffv_dens);

  ffv_ffv_vec_test(vecffv_prob, vecffv_dens);
  ffv_ffv_vec_test(vecffv_prob, row_vecffv_dens);
  ffv_ffv_vec_test(vecffv_prob, std_vecffv_dens);

  ffv_ffv_test(row_vecffv_prob, vecffv_dens);
  ffv_ffv_test(row_vecffv_prob, row_vecffv_dens);
  ffv_ffv_test(row_vecffv_prob, std_vecffv_dens);

  ffv_ffv_test(std_vecffv_prob, vecffv_dens);
  ffv_ffv_test(std_vecffv_prob, row_vecffv_dens);
  ffv_ffv_test(std_vecffv_prob, std_vecffv_dens);
}

template <typename T_a, typename T_b>
void ffv_d_test(T_a a, T_b b) {
  a[0].val_ = 0.514;
  a[1].val_ = 0.284;
  a[2].val_ = 0.112;
  a[3].val_ = 0.090;
  a[0].d_ = 1.0;
  a[1].d_ = 1.0;
  a[2].d_ = 1.0;
  a[3].d_ = 1.0;
  a[0].val_.d_ = 1.0;
  a[1].val_.d_ = 1.0;
  a[2].val_.d_ = 1.0;
  a[3].val_.d_ = 1.0;

  b[0] = -3.581;
  b[1] = -8.114;
  b[2] = -11.215;
  b[3] = -5.658;

  fvar<fvar<var> > out = log_mix(a, b);

  EXPECT_FLOAT_EQ(out.val_.val_.val(), -4.218931574);
  EXPECT_FLOAT_EQ(out.d_.val_.val(), 2.150968236);
  out.d_.val_.grad();

  vector_d prob_deriv(4);
  prob_deriv << 1.892562202198, 0.020341982471, 0.000915474153, 0.237148577149;

  for (int i = 0; i < 4; ++i) {
    EXPECT_FLOAT_EQ(a[i].d_.val_.adj(), prob_deriv[i]);
  }
}

template <typename T_a, typename T_b>
void ffv_d_vec_test(T_a a, T_b b1) {
  a[0].val_ = 0.514;
  a[1].val_ = 0.284;
  a[2].val_ = 0.112;
  a[3].val_ = 0.090;
  a[0].d_ = 1.0;
  a[1].d_ = 1.0;
  a[2].d_ = 1.0;
  a[3].d_ = 1.0;
  a[0].val_.d_ = 1.0;
  a[1].val_.d_ = 1.0;
  a[2].val_.d_ = 1.0;
  a[3].val_.d_ = 1.0;

  b1[0] = -3.581;
  b1[1] = -8.114;
  b1[2] = -11.215;
  b1[3] = -5.658;

  T_b b2(4), b3(4);

  b2[0] = -8.594;
  b2[1] = -3.251;
  b2[2] = -7.281;
  b2[3] = -3.556;

  b3[0] = -11.554;
  b3[1] = -6.628;
  b3[2] = -15.229;
  b3[3] = -9.561;

  std::vector<T_b> c{b1, b2, b3};

  fvar<fvar<var> > out = log_mix(a, c);

  EXPECT_FLOAT_EQ(out.val_.val_.val(), -16.36331174);
  EXPECT_FLOAT_EQ(out.d_.val_.val(), 10.73648479);
  out.d_.val_.grad();

  vector_d prob_deriv(4);
  prob_deriv << 1.930840738468, 6.257235890773, 0.051642416014, 2.496765742829;

  for (int i = 0; i < 4; ++i) {
    EXPECT_FLOAT_EQ(a[i].d_.val_.adj(), prob_deriv[i]);
  }
}

TEST(AgradMixMatrixLogMix, ffv_d) {
  vector_ffv vecffv_prob(4);
  vector_d vecd_dens(4);
  row_vector_ffv row_vecffv_prob(4);
  row_vector_d row_vecd_dens(4);
  std::vector<fvar<fvar<var> > > std_vecffv_prob(4);
  std::vector<double> std_vecd_dens(4);

  ffv_d_test(vecffv_prob, vecd_dens);
  ffv_d_test(vecffv_prob, row_vecd_dens);
  ffv_d_test(vecffv_prob, std_vecd_dens);

  ffv_d_vec_test(vecffv_prob, vecd_dens);
  ffv_d_vec_test(vecffv_prob, row_vecd_dens);
  ffv_d_vec_test(vecffv_prob, std_vecd_dens);

  ffv_d_test(row_vecffv_prob, vecd_dens);
  ffv_d_test(row_vecffv_prob, row_vecd_dens);
  ffv_d_test(row_vecffv_prob, std_vecd_dens);

  ffv_d_vec_test(row_vecffv_prob, vecd_dens);
  ffv_d_vec_test(row_vecffv_prob, row_vecd_dens);
  ffv_d_vec_test(row_vecffv_prob, std_vecd_dens);

  ffv_d_test(std_vecffv_prob, vecd_dens);
  ffv_d_test(std_vecffv_prob, row_vecd_dens);
  ffv_d_test(std_vecffv_prob, std_vecd_dens);

  ffv_d_vec_test(std_vecffv_prob, vecd_dens);
  ffv_d_vec_test(std_vecffv_prob, row_vecd_dens);
  ffv_d_vec_test(std_vecffv_prob, std_vecd_dens);
}

template <typename T_a, typename T_b>
void d_ffv_test(T_a a, T_b b) {
  a[0] = 0.514;
  a[1] = 0.284;
  a[2] = 0.112;
  a[3] = 0.090;

  b[0].val_ = -3.581;
  b[1].val_ = -8.114;
  b[2].val_ = -11.215;
  b[3].val_ = -5.658;
  b[0].d_ = 1.0;
  b[1].d_ = 1.0;
  b[2].d_ = 1.0;
  b[3].d_ = 1.0;
  b[0].val_.d_ = 1.0;
  b[1].val_.d_ = 1.0;
  b[2].val_.d_ = 1.0;
  b[3].val_.d_ = 1.0;

  fvar<fvar<var> > out = log_mix(a, b);

  EXPECT_FLOAT_EQ(out.val_.val_.val(), -4.218931574);
  EXPECT_FLOAT_EQ(out.d_.val_.val(), 1.0);
  out.d_.val_.grad();

  vector_d dens_deriv(4);
  dens_deriv << 0.97277697193, 0.00577712302, 0.00010253311, 0.02134337194;

  for (int i = 0; i < 4; ++i) {
    EXPECT_FLOAT_EQ(b[i].d_.val_.adj(), dens_deriv[i]);
  }
}

template <typename T_a, typename T_b>
void d_ffv_vec_test(T_a a, T_b b1) {
  a[0] = 0.514;
  a[1] = 0.284;
  a[2] = 0.112;
  a[3] = 0.090;

  b1[0].val_ = -3.581;
  b1[1].val_ = -8.114;
  b1[2].val_ = -11.215;
  b1[3].val_ = -5.658;

  b1[0].d_ = 1.0;
  b1[1].d_ = 1.0;
  b1[2].d_ = 1.0;
  b1[3].d_ = 1.0;
  b1[0].val_.d_ = 1.0;
  b1[1].val_.d_ = 1.0;
  b1[2].val_.d_ = 1.0;
  b1[3].val_.d_ = 1.0;

  T_b b2(4), b3(4);

  b2[0].val_ = -8.594;
  b2[1].val_ = -3.251;
  b2[2].val_ = -7.281;
  b2[3].val_ = -3.556;

  b2[0].d_ = 1.0;
  b2[1].d_ = 1.0;
  b2[2].d_ = 1.0;
  b2[3].d_ = 1.0;
  b2[0].val_.d_ = 1.0;
  b2[1].val_.d_ = 1.0;
  b2[2].val_.d_ = 1.0;
  b2[3].val_.d_ = 1.0;

  b3[0].val_ = -11.554;
  b3[1].val_ = -6.628;
  b3[2].val_ = -15.229;
  b3[3].val_ = -9.561;

  b3[0].d_ = 1.0;
  b3[1].d_ = 1.0;
  b3[2].d_ = 1.0;
  b3[3].d_ = 1.0;
  b3[0].val_.d_ = 1.0;
  b3[1].val_.d_ = 1.0;
  b3[2].val_.d_ = 1.0;
  b3[3].val_.d_ = 1.0;

  std::vector<T_b> c{b1, b2, b3};

  fvar<fvar<var> > out = log_mix(a, c);

  EXPECT_FLOAT_EQ(out.val_.val_.val(), -16.36331174);
  EXPECT_FLOAT_EQ(out.d_.val_.val(), 3.0);
  out.d_.val_.grad();

  vector_d dens1_deriv(4);
  dens1_deriv << 0.97277697193, 0.00577712302, 0.00010253311, 0.02134337194;
  vector_d dens2_deriv(4);
  dens2_deriv << 0.00692718708, 0.80047459791, 0.00561100267, 0.18698721234;
  vector_d dens3_deriv(4);
  dens3_deriv << 0.01274798056, 0.97080327205, 0.00007041482, 0.01637833257;

  for (int i = 0; i < 4; ++i) {
    EXPECT_FLOAT_EQ(c[0][i].d_.val_.adj(), dens1_deriv[i]);
    EXPECT_FLOAT_EQ(c[1][i].d_.val_.adj(), dens2_deriv[i]);
    EXPECT_FLOAT_EQ(c[2][i].d_.val_.adj(), dens3_deriv[i]);
  }
}

TEST(AgradMixMatrixLogMix, d_ffv) {
  vector_d vecd_prob(4);
  vector_ffv vecffv_dens(4);
  row_vector_d row_vecd_prob(4);
  row_vector_ffv row_vecffv_dens(4);
  std::vector<double> std_vecd_prob(4);
  std::vector<fvar<fvar<var> > > std_vecffv_dens(4);

  d_ffv_test(vecd_prob, vecffv_dens);
  d_ffv_test(vecd_prob, row_vecffv_dens);
  d_ffv_test(vecd_prob, std_vecffv_dens);

  d_ffv_vec_test(vecd_prob, vecffv_dens);
  d_ffv_vec_test(vecd_prob, row_vecffv_dens);
  d_ffv_vec_test(vecd_prob, std_vecffv_dens);

  d_ffv_test(row_vecd_prob, vecffv_dens);
  d_ffv_test(row_vecd_prob, row_vecffv_dens);
  d_ffv_test(row_vecd_prob, std_vecffv_dens);

  d_ffv_vec_test(row_vecd_prob, vecffv_dens);
  d_ffv_vec_test(row_vecd_prob, row_vecffv_dens);
  d_ffv_vec_test(row_vecd_prob, std_vecffv_dens);

  d_ffv_test(std_vecd_prob, vecffv_dens);
  d_ffv_test(std_vecd_prob, row_vecffv_dens);
  d_ffv_test(std_vecd_prob, std_vecffv_dens);

  d_ffv_vec_test(std_vecd_prob, vecffv_dens);
  d_ffv_vec_test(std_vecd_prob, row_vecffv_dens);
  d_ffv_vec_test(std_vecd_prob, std_vecffv_dens);
}

TEST(AgradMixMatrixLogMix, fv_fv_old) {
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

TEST(AgradMixMatrixLogMix, fv_d_old) {
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

TEST(AgradMixMatrixLogMix, d_fv_old) {
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

TEST(AgradMixMatrixLogMix, ffv_ffv_old) {
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

TEST(AgradMixMatrixLogMix, ffv_d_old) {
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

TEST(AgradMixMatrixLogMix, d_ffv_old) {
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
