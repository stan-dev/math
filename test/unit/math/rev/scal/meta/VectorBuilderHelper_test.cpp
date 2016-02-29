#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, VectorBuilderHelper_false_true) {
  using stan::VectorBuilderHelper;
  using stan::math::var;
  using stan::length;

  var a_var(1);

  VectorBuilderHelper<double,false,true> dvv1(length(a_var));
  EXPECT_THROW(dvv1[0], std::logic_error);
}

TEST(MetaTraits, VectorBuilderHelper_true_true) {
  using stan::VectorBuilderHelper;
  using stan::math::var;
  using stan::length;

  var a_var(1);

  VectorBuilderHelper<double,true,true> dvv1(length(a_var));
  dvv1[0] = 0.0;
  EXPECT_FLOAT_EQ(0.0, dvv1[0]);
}
