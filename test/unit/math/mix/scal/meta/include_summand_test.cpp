#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::math::fvar;
using stan::math::include_summand;
using stan::math::var;

TEST(MetaTraits, IncludeSummandProptoTrueFvarVar) {
  EXPECT_TRUE((include_summand<true, fvar<var> >::value));
}

TEST(MetaTraits, IncludeSummandProptoTrueFvarFvarVar) {
  EXPECT_TRUE((include_summand<true, fvar<fvar<var> > >::value));
}

TEST(MetaTraits, IncludeSummandProtoTrueFvarVarTen) {
  EXPECT_TRUE((include_summand<true, double, fvar<var>, int, fvar<var>, double,
                               double, int, int, fvar<var>, int>::value));
}

TEST(MetaTraits, IncludeSummandProtoTrueFvarFvarVarTen) {
  EXPECT_TRUE((
      include_summand<true, double, fvar<fvar<var> >, int, fvar<fvar<var> >,
                      double, double, int, int, fvar<fvar<var> >, int>::value));
}
