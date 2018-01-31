#include <stan/math/rev/mat.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(AgradRevMatrix, log_mix_avec) {
  AVEC prob{0.15, 0.20, 0.40, 0.25};
  AVEC dens{-2.15, -3.89, -2.18, -8.82};

  AVAR out = stan::math::log_mix(prob, dens);

  out.grad();

  EXPECT_FLOAT_EQ(out.val(), -2.70582405);

  EXPECT_FLOAT_EQ(prob[0].adj(), 1.7433770185);
  EXPECT_FLOAT_EQ(prob[1].adj(), 0.3059982327);
  EXPECT_FLOAT_EQ(prob[2].adj(), 1.6918524409);
  EXPECT_FLOAT_EQ(prob[3].adj(), 0.0022112972);

  EXPECT_FLOAT_EQ(dens[0].adj(), 0.2615065527);
  EXPECT_FLOAT_EQ(dens[1].adj(), 0.0611996465);
  EXPECT_FLOAT_EQ(dens[2].adj(), 0.6767409763);
  EXPECT_FLOAT_EQ(dens[3].adj(), 0.0005528243);
}

TEST(AgradRevMatrix, log_mix_vector_v) {
  using stan::math::vector_v;

  vector_v prob(4), dens(4);
  prob << 0.13, 0.22, 0.38, 0.27;
  dens << -3.15, -0.21, -10.55, -7.24;

  AVAR out = stan::math::log_mix(prob, dens);

  out.grad();

  EXPECT_FLOAT_EQ(out.val(), -1.69226023);

  EXPECT_FLOAT_EQ(prob[0].adj(), 0.2327617758);
  EXPECT_FLOAT_EQ(prob[1].adj(), 4.4028859801);
  EXPECT_FLOAT_EQ(prob[2].adj(), 0.0001422763);
  EXPECT_FLOAT_EQ(prob[3].adj(), 0.0038962537);

  EXPECT_FLOAT_EQ(dens[0].adj(), 0.0302590308);
  EXPECT_FLOAT_EQ(dens[1].adj(), 0.9686349156);
  EXPECT_FLOAT_EQ(dens[2].adj(), 5.406498E-05);
  EXPECT_FLOAT_EQ(dens[3].adj(), 0.0010519885);
}

TEST(AgradRevMatrix, log_mix_row_vector_v) {
  using stan::math::row_vector_v;

  row_vector_v prob(4), dens(4);
  prob << 0.03, 0.21, 0.63, 0.13;
  dens << -19.41, -8.14, -2.18, -9.13;

  AVAR out = stan::math::log_mix(prob, dens);

  out.grad();

  EXPECT_FLOAT_EQ(out.val(), -2.64097823);

  EXPECT_FLOAT_EQ(prob[0].adj(), 5.215625E-08);
  EXPECT_FLOAT_EQ(prob[1].adj(), 0.0040907712);
  EXPECT_FLOAT_EQ(prob[2].adj(), 1.5856243363);
  EXPECT_FLOAT_EQ(prob[3].adj(), 0.0015200352);

  EXPECT_FLOAT_EQ(dens[0].adj(), 1.564688E-09);
  EXPECT_FLOAT_EQ(dens[1].adj(), 0.0008590619);
  EXPECT_FLOAT_EQ(dens[2].adj(), 0.9989433319);
  EXPECT_FLOAT_EQ(dens[3].adj(), 0.0001976046);
}

TEST(AgradRevMatrix, log_mix_std_vector_v) {
  using stan::math::var;

  std::vector<var> prob(4);
  prob[0] = 0.03;
  prob[1] = 0.21;
  prob[2] = 0.63;
  prob[3] = 0.13;
  std::vector<var> dens(4);
  dens[0] = -19.41;
  dens[1] = -8.14;
  dens[2] = -2.18;
  dens[3] = -9.13;

  AVAR out = stan::math::log_mix(prob, dens);

  out.grad();

  EXPECT_FLOAT_EQ(out.val(), -2.64097823);

  EXPECT_FLOAT_EQ(prob[0].adj(), 5.215625E-08);
  EXPECT_FLOAT_EQ(prob[1].adj(), 0.0040907712);
  EXPECT_FLOAT_EQ(prob[2].adj(), 1.5856243363);
  EXPECT_FLOAT_EQ(prob[3].adj(), 0.0015200352);

  EXPECT_FLOAT_EQ(dens[0].adj(), 1.564688E-09);
  EXPECT_FLOAT_EQ(dens[1].adj(), 0.0008590619);
  EXPECT_FLOAT_EQ(dens[2].adj(), 0.9989433319);
  EXPECT_FLOAT_EQ(dens[3].adj(), 0.0001976046);
}
