#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/util.hpp>

TEST(ProbInternalMath, gradF32_fd1) {
  using stan::math::fvar;

  fvar<double> a = 1.0;
  a.d_ = 1.0;
  fvar<double> b = 31.0;
  fvar<double> c = -27.0;
  fvar<double> d = 19.0;
  fvar<double> e = -41.0;
  fvar<double> z = 1.0;
  fvar<double> g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_,1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_,1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_,1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_,1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_,1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_,1e-5);

  EXPECT_NEAR(41.01553475870347475023037358640582917147051389292162474016745,
              g[0].d_,1e-5);
}

TEST(ProbInternalMath, gradF32_fd2) {
  using stan::math::fvar;

  fvar<double> a = 1.0;
  fvar<double> b = 31.0;
  b.d_ = 1.0;
  fvar<double> c = -27.0;
  fvar<double> d = 19.0;
  fvar<double> e = -41.0;
  fvar<double> z = 1.0;
  fvar<double> g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_,1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_,1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_,1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_,1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_,1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_,1e-5);

  EXPECT_NEAR(0.342454543339724329115552426438001592723143365030924900588111,
              g[1].d_,1e-5);
}

TEST(ProbInternalMath, gradF32_fd3) {
  using stan::math::fvar;

  fvar<double> a = 1.0;
  fvar<double> b = 31.0;
  fvar<double> c = -27.0;
  c.d_ = 1.0;
  fvar<double> d = 19.0;
  fvar<double> e = -41.0;
  fvar<double> z = 1.0;
  fvar<double> g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_,1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_,1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_,1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_,1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_,1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_,1e-5);

  EXPECT_NEAR(0.90986472078762437,
              g[2].d_,1e-5);
}
TEST(ProbInternalMath, gradF32_fd4) {
  using stan::math::fvar;

  fvar<double> a = 1.0;
  fvar<double> b = 31.0;
  fvar<double> c = -27.0;
  fvar<double> d = 19.0;
  d.d_ = 1.0;
  fvar<double> e = -41.0;
  fvar<double> z = 1.0;
  fvar<double> g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_,1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_,1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_,1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_,1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_,1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_,1e-5);

  EXPECT_NEAR(1.047024959065504556655904003595645684444382378830047020988218,
              g[3].d_,1e-5);
}
TEST(ProbInternalMath, gradF32_fd5) {
  using stan::math::fvar;

  fvar<double> a = 1.0;
  fvar<double> b = 31.0;
  fvar<double> c = -27.0;
  fvar<double> d = 19.0;
  fvar<double> e = -41.0;
  e.d_ = 1.0;
  fvar<double> z = 1.0;
  fvar<double> g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_,1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_,1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_,1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_,1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_,1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_,1e-5);

  EXPECT_NEAR(0.415359887777218792995404669803015764396172842233556866773418,
              g[4].d_,1e-5);
}
TEST(ProbInternalMath, gradF32_fd6) {
  using stan::math::fvar;

  fvar<double> a = 1.0;
  fvar<double> b = 31.0;
  fvar<double> c = -27.0;
  fvar<double> d = 19.0;
  fvar<double> e = -41.0;
  fvar<double> z = 1.0;
  z.d_ = 1.0;
  fvar<double> g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_,1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_,1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_,1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_,1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_,1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_,1e-5);

  EXPECT_NEAR(424.5724606148232594702100102534498155985480235827583548085963,
              g[5].d_,1e-5);
}

TEST(ProbInternalMath, gradF32_ffd_2ndderiv1) {
  using stan::math::fvar;

  fvar<fvar<double> > a = 1.0;
  a.d_ = 1.0;
  fvar<fvar<double> > b = 31.0;
  fvar<fvar<double> > c = -27.0;
  fvar<fvar<double> > d = 19.0;
  fvar<fvar<double> > e = -41.0;
  fvar<fvar<double> > z = 1.0;
  fvar<fvar<double> > g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val_,1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val_,1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val_,1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val_,1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val_,1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val_,1e-5);

  EXPECT_NEAR(41.01553475870347475023037358640582917147051389292162474016745,
              g[0].d_.val_,1e-5);
}

TEST(ProbInternalMath, gradF32_ffd_2ndderiv2) {
  using stan::math::fvar;

  fvar<fvar<double> > a = 1.0;
  fvar<fvar<double> > b = 31.0;
  b.d_ = 1.0;
  fvar<fvar<double> > c = -27.0;
  fvar<fvar<double> > d = 19.0;
  fvar<fvar<double> > e = -41.0;
  fvar<fvar<double> > z = 1.0;
  fvar<fvar<double> > g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val_,1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val_,1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val_,1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val_,1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val_,1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val_,1e-5);

  EXPECT_NEAR(0.342454543339724329115552426438001592723143365030924900588111,
              g[1].d_.val_,1e-5);
}

TEST(ProbInternalMath, gradF32_ffd_2ndderiv3) {
  using stan::math::fvar;

  fvar<fvar<double> > a = 1.0;
  fvar<fvar<double> > b = 31.0;
  fvar<fvar<double> > c = -27.0;
  c.d_ = 1.0;
  fvar<fvar<double> > d = 19.0;
  fvar<fvar<double> > e = -41.0;
  fvar<fvar<double> > z = 1.0;
  fvar<fvar<double> > g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val_,1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val_,1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val_,1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val_,1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val_,1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val_,1e-5);

  EXPECT_NEAR(0.90986472078762437,
              g[2].d_.val_,1e-5);
}
TEST(ProbInternalMath, gradF32_ffd_2ndderiv4) {
  using stan::math::fvar;

  fvar<fvar<double> > a = 1.0;
  fvar<fvar<double> > b = 31.0;
  fvar<fvar<double> > c = -27.0;
  fvar<fvar<double> > d = 19.0;
  d.d_ = 1.0;
  fvar<fvar<double> > e = -41.0;
  fvar<fvar<double> > z = 1.0;
  fvar<fvar<double> > g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val_,1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val_,1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val_,1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val_,1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val_,1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val_,1e-5);

  EXPECT_NEAR(1.047024959065504556655904003595645684444382378830047020988218,
              g[3].d_.val_,1e-5);
}
TEST(ProbInternalMath, gradF32_ffd_2ndderiv5) {
  using stan::math::fvar;

  fvar<fvar<double> > a = 1.0;
  fvar<fvar<double> > b = 31.0;
  fvar<fvar<double> > c = -27.0;
  fvar<fvar<double> > d = 19.0;
  fvar<fvar<double> > e = -41.0;
  e.d_ = 1.0;
  fvar<fvar<double> > z = 1.0;
  fvar<fvar<double> > g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val_,1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val_,1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val_,1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val_,1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val_,1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val_,1e-5);

  EXPECT_NEAR(0.415359887777218792995404669803015764396172842233556866773418,
              g[4].d_.val_,1e-5);
}
TEST(ProbInternalMath, gradF32_ffd_2ndderiv6) {
  using stan::math::fvar;

  fvar<fvar<double> > a = 1.0;
  fvar<fvar<double> > b = 31.0;
  fvar<fvar<double> > c = -27.0;
  fvar<fvar<double> > d = 19.0;
  fvar<fvar<double> > e = -41.0;
  fvar<fvar<double> > z = 1.0;
  z.d_ = 1.0;
  fvar<fvar<double> > g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val_,1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val_,1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val_,1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val_,1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val_,1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val_,1e-5);

  EXPECT_NEAR(424.5724606148232594702100102534498155985480235827583548085963,
              g[5].d_.val_,1e-5);
}

TEST(ProbInternalMath, gradF32_ffd_3rdderiv1) {
  using stan::math::fvar;

  fvar<fvar<double> > a = 1.0;
  a.d_ = 1.0;
  a.val_.d_ = 1.0;
  fvar<fvar<double> > b = 31.0;
  fvar<fvar<double> > c = -27.0;
  fvar<fvar<double> > d = 19.0;
  fvar<fvar<double> > e = -41.0;
  fvar<fvar<double> > z = 1.0;
  fvar<fvar<double> > g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val_,1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val_,1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val_,1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val_,1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val_,1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val_,1e-5);

  EXPECT_NEAR(65.599396543196708101135082478886528455,g[0].d_.d_,1e-5);
}

TEST(ProbInternalMath, gradF32_ffd_3rdderiv2) {
  using stan::math::fvar;

  fvar<fvar<double> > a = 1.0;
  fvar<fvar<double> > b = 31.0;
  b.d_ = 1.0;
  b.val_.d_ = 1.0;
  fvar<fvar<double> > c = -27.0;
  fvar<fvar<double> > d = 19.0;
  fvar<fvar<double> > e = -41.0;
  fvar<fvar<double> > z = 1.0;
  fvar<fvar<double> > g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val_,1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val_,1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val_,1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val_,1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val_,1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val_,1e-5);

  EXPECT_NEAR(0.074524365399251673905543999863038234,g[1].d_.d_,1e-5);
}

TEST(ProbInternalMath, gradF32_ffd_3rdderiv3) {
  using stan::math::fvar;

  fvar<fvar<double> > a = 1.0;
  fvar<fvar<double> > b = 31.0;
  fvar<fvar<double> > c = -27.0;
  c.d_ = 1.0;
  c.val_.d_ = 1.0;
  fvar<fvar<double> > d = 19.0;
  fvar<fvar<double> > e = -41.0;
  fvar<fvar<double> > z = 1.0;
  fvar<fvar<double> > g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val_,1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val_,1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val_,1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val_,1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val_,1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val_,1e-5);

  EXPECT_NEAR(-0.4025421605307411,g[2].d_.d_,1e-5);
}
TEST(ProbInternalMath, gradF32_ffd_3rdderiv4) {
  using stan::math::fvar;

  fvar<fvar<double> > a = 1.0;
  fvar<fvar<double> > b = 31.0;
  fvar<fvar<double> > c = -27.0;
  fvar<fvar<double> > d = 19.0;
  d.d_ = 1.0;
  d.val_.d_ = 1.0;
  fvar<fvar<double> > e = -41.0;
  fvar<fvar<double> > z = 1.0;
  fvar<fvar<double> > g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val_,1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val_,1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val_,1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val_,1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val_,1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val_,1e-5);

  EXPECT_NEAR(-0.505769456958641747831864908555691738,g[3].d_.d_,1e-5);
}
TEST(ProbInternalMath, gradF32_ffd_3rdderiv5) {
  using stan::math::fvar;

  fvar<fvar<double> > a = 1.0;
  fvar<fvar<double> > b = 31.0;
  fvar<fvar<double> > c = -27.0;
  fvar<fvar<double> > d = 19.0;
  fvar<fvar<double> > e = -41.0;
  e.d_ = 1.0;
  e.val_.d_ = 1.0;
  fvar<fvar<double> > z = 1.0;
  fvar<fvar<double> > g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val_,1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val_,1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val_,1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val_,1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val_,1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val_,1e-5);

  EXPECT_NEAR(0.143334459434983770772868891143925349,g[4].d_.d_,1e-5);
}
TEST(ProbInternalMath, gradF32_ffd_3rdderiv6) {
  using stan::math::fvar;

  fvar<fvar<double> > a = 1.0;
  fvar<fvar<double> > b = 31.0;
  fvar<fvar<double> > c = -27.0;
  fvar<fvar<double> > d = 19.0;
  fvar<fvar<double> > e = -41.0;
  fvar<fvar<double> > z = 1.0;
  z.d_ = 1.0;
  z.val_.d_ = 1.0;
  fvar<fvar<double> > g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val_,1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val_,1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val_,1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val_,1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val_,1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val_,1e-5);

  EXPECT_NEAR(3464.701930495754952696665090654133564817,g[5].d_.d_,1e-5);
}

TEST(ProbInternalMath, gradF32_fv_1stderiv1) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 1.0;
  a.d_ = 1.0;
  fvar<var> b = 31.0;
  fvar<var> c = -27.0;
  fvar<var> d = 19.0;
  fvar<var> e = -41.0;
  fvar<var> z = 1.0;
  fvar<var> g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val(),1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val(),1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val(),1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val(),1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val(),1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val(),1e-5);

  EXPECT_NEAR(41.01553475870347475023037358640582917147051389292162474016745,
              g[0].d_.val(),1e-5);

  AVEC y1 = createAVEC(a.val_);
  VEC grad1;
  g[0].val_.grad(y1,grad1);
  EXPECT_NEAR(41.01553475870347475023037358640582917147051389292162474016745,grad1[0],1e-5);
}

TEST(ProbInternalMath, gradF32_fv_1stderiv2) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 1.0;
  fvar<var> b = 31.0;
  b.d_ = 1.0;
  fvar<var> c = -27.0;
  fvar<var> d = 19.0;
  fvar<var> e = -41.0;
  fvar<var> z = 1.0;
  fvar<var> g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val(),1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val(),1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val(),1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val(),1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val(),1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val(),1e-5);

  EXPECT_NEAR(0.342454543339724329115552426438001592723143365030924900588111,
              g[1].d_.val(),1e-5);

  AVEC y1 = createAVEC(b.val_);
  VEC grad1;
  g[1].val_.grad(y1,grad1);
  EXPECT_NEAR(0.342454543339724329115552426438001592723143365030924900588111,grad1[0],1e-5);
}

TEST(ProbInternalMath, gradF32_fv_1stderiv3) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 1.0;
  fvar<var> b = 31.0;
  fvar<var> c = -27.0;
  c.d_ = 1.0;
  fvar<var> d = 19.0;
  fvar<var> e = -41.0;
  fvar<var> z = 1.0;
  fvar<var> g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val(),1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val(),1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val(),1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val(),1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val(),1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val(),1e-5);

  EXPECT_NEAR(0.90986472078762437,
              g[2].d_.val(),1e-5);

  AVEC y1 = createAVEC(c.val_);
  VEC grad1;
  g[2].val_.grad(y1,grad1);
  EXPECT_NEAR(0.90986472078762437,grad1[0],1e-5);
}
TEST(ProbInternalMath, gradF32_fv_1stderiv4) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 1.0;
  fvar<var> b = 31.0;
  fvar<var> c = -27.0;
  fvar<var> d = 19.0;
  d.d_ = 1.0;
  fvar<var> e = -41.0;
  fvar<var> z = 1.0;
  fvar<var> g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val(),1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val(),1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val(),1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val(),1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val(),1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val(),1e-5);

  EXPECT_NEAR(1.047024959065504556655904003595645684444382378830047020988218,
              g[3].d_.val(),1e-5);

  AVEC y1 = createAVEC(d.val_);
  VEC grad1;
  g[3].val_.grad(y1,grad1);
  EXPECT_NEAR(1.047024959065504556655904003595645684444382378830047020988218,grad1[0],1e-5);
}
TEST(ProbInternalMath, gradF32_fv_1stderiv5) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 1.0;
  fvar<var> b = 31.0;
  fvar<var> c = -27.0;
  fvar<var> d = 19.0;
  fvar<var> e = -41.0;
  e.d_ = 1.0;
  fvar<var> z = 1.0;
  fvar<var> g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val(),1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val(),1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val(),1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val(),1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val(),1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val(),1e-5);

  EXPECT_NEAR(0.415359887777218792995404669803015764396172842233556866773418,
              g[4].d_.val(),1e-5);

  AVEC y1 = createAVEC(e.val_);
  VEC grad1;
  g[4].val_.grad(y1,grad1);
  EXPECT_NEAR(0.415359887777218792995404669803015764396172842233556866773418,grad1[0],1e-5);
}
TEST(ProbInternalMath, gradF32_fv_1stderiv6) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 1.0;
  fvar<var> b = 31.0;
  fvar<var> c = -27.0;
  fvar<var> d = 19.0;
  fvar<var> e = -41.0;
  fvar<var> z = 1.0;
  z.d_ = 1.0;
  fvar<var> g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val(),1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val(),1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val(),1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val(),1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val(),1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val(),1e-5);

  EXPECT_NEAR(424.5724606148232594702100102534498155985480235827583548085963,
              g[5].d_.val(),1e-5);

  AVEC y1 = createAVEC(z.val_);
  VEC grad1;
  g[5].val_.grad(y1,grad1);
  EXPECT_NEAR(424.5724606148232594702100102534498155985480235827583548085963,grad1[0],1e-5);
}

TEST(ProbInternalMath, gradF32_fv_2ndderiv1) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 1.0;
  a.d_ = 1.0;
  fvar<var> b = 31.0;
  fvar<var> c = -27.0;
  fvar<var> d = 19.0;
  fvar<var> e = -41.0;
  fvar<var> z = 1.0;
  fvar<var> g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val(),1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val(),1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val(),1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val(),1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val(),1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val(),1e-5);

  EXPECT_NEAR(41.01553475870347475023037358640582917147051389292162474016745,
              g[0].d_.val(),1e-5);

  AVEC y1 = createAVEC(a.val_);
  VEC grad1;
  g[0].d_.grad(y1,grad1);
  EXPECT_NEAR(65.599396543196708101135082478886528455,grad1[0],1e-5);
}

TEST(ProbInternalMath, gradF32_fv_2ndderiv2) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 1.0;
  fvar<var> b = 31.0;
  b.d_ = 1.0;
  fvar<var> c = -27.0;
  fvar<var> d = 19.0;
  fvar<var> e = -41.0;
  fvar<var> z = 1.0;
  fvar<var> g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val(),1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val(),1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val(),1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val(),1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val(),1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val(),1e-5);

  EXPECT_NEAR(0.342454543339724329115552426438001592723143365030924900588111,
              g[1].d_.val(),1e-5);

  AVEC y1 = createAVEC(b.val_);
  VEC grad1;
  g[1].d_.grad(y1,grad1);
  EXPECT_NEAR(0.074524365399251673905543999863038234,grad1[0],1e-5);
}

TEST(ProbInternalMath, gradF32_fv_2ndderiv3) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 1.0;
  fvar<var> b = 31.0;
  fvar<var> c = -27.0;
  c.d_ = 1.0;
  fvar<var> d = 19.0;
  fvar<var> e = -41.0;
  fvar<var> z = 1.0;
  fvar<var> g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val(),1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val(),1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val(),1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val(),1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val(),1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val(),1e-5);

  EXPECT_NEAR(0.90986472078762437,
              g[2].d_.val(),1e-5);

  AVEC y1 = createAVEC(c.val_);
  VEC grad1;
  g[2].d_.grad(y1,grad1);
  EXPECT_NEAR(-0.4025421605307411,grad1[0],1e-5);
}
TEST(ProbInternalMath, gradF32_fv_2ndderiv4) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 1.0;
  fvar<var> b = 31.0;
  fvar<var> c = -27.0;
  fvar<var> d = 19.0;
  d.d_ = 1.0;
  fvar<var> e = -41.0;
  fvar<var> z = 1.0;
  fvar<var> g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val(),1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val(),1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val(),1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val(),1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val(),1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val(),1e-5);

  EXPECT_NEAR(1.047024959065504556655904003595645684444382378830047020988218,
              g[3].d_.val(),1e-5);

  AVEC y1 = createAVEC(d.val_);
  VEC grad1;
  g[3].d_.grad(y1,grad1);
  EXPECT_NEAR(-0.505769456958641747831864908555691738,grad1[0],1e-5);
}
TEST(ProbInternalMath, gradF32_fv_2ndderiv5) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 1.0;
  fvar<var> b = 31.0;
  fvar<var> c = -27.0;
  fvar<var> d = 19.0;
  fvar<var> e = -41.0;
  e.d_ = 1.0;
  fvar<var> z = 1.0;
  fvar<var> g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val(),1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val(),1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val(),1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val(),1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val(),1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val(),1e-5);

  EXPECT_NEAR(0.415359887777218792995404669803015764396172842233556866773418,
              g[4].d_.val(),1e-5);

  AVEC y1 = createAVEC(e.val_);
  VEC grad1;
  g[4].d_.grad(y1,grad1);
  EXPECT_NEAR(0.143334459434983770772868891143925349,grad1[0],1e-5);
}
TEST(ProbInternalMath, gradF32_fv_2ndderiv6) {
  using stan::math::fvar;
  using stan::math::var;

  fvar<var> a = 1.0;
  fvar<var> b = 31.0;
  fvar<var> c = -27.0;
  fvar<var> d = 19.0;
  fvar<var> e = -41.0;
  fvar<var> z = 1.0;
  z.d_ = 1.0;
  fvar<var> g[6];

  stan::math::grad_F32(g,a,b,c,d,e,z);

  EXPECT_NEAR(22.95829816018250585416491584581112223816561212219172212450836,
              g[0].val_.val(),1e-5);
  EXPECT_NEAR(1.740056451478897241488082512854205170874142224663970334770766,
              g[1].val_.val(),1e-5);
  EXPECT_NEAR(-2.6052400424887519,
              g[2].val_.val(),1e-5);
  EXPECT_NEAR(-2.69297893625022464634137707353872105148995696636392529052847,
              g[3].val_.val(),1e-5);
  EXPECT_NEAR(1.606519030743225019685406202522547937181609777973133379643512,
              g[4].val_.val(),1e-5);
  EXPECT_NEAR(59.65791128638963495870649759618269069328482053152963195380560,
              g[5].val_.val(),1e-5);

  EXPECT_NEAR(424.5724606148232594702100102534498155985480235827583548085963,
              g[5].d_.val(),1e-5);

  AVEC y1 = createAVEC(z.val_);
  VEC grad1;
  g[5].d_.grad(y1,grad1);
  EXPECT_NEAR(3464.701930495754952696665090654133564817,grad1[0],1e-5);
}
