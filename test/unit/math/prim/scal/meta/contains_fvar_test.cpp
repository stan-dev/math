#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <boost/type_traits.hpp>

TEST(MetaTraits,containsFvar) {
  using stan::contains_fvar;
  EXPECT_FALSE(contains_fvar<double>::value);
}
