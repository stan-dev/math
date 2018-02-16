#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, VectorBuilderHelper_false_true) {
  using stan::VectorBuilderHelper;
  using stan::length;
  using stan::math::var;

  var a_var(1);

  VectorBuilderHelper<double, false, true> dvv1(length(a_var));
  EXPECT_THROW(dvv1[0], std::logic_error);
  EXPECT_THROW(dvv1.data(), std::logic_error);
}

TEST(MetaTraits, VectorBuilderHelper_true_true) {
  using stan::VectorBuilderHelper;
  using stan::length;
  using stan::math::var;

  var a_var(1);

  VectorBuilderHelper<double, true, true> dvv1(length(a_var));
  EXPECT_THROW(dvv1[0], std::logic_error)
      << "This uses the default template; if the arr version is included, "
      << "it will use the template specialization.";
}
