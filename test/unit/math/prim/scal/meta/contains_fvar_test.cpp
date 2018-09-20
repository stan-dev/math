#include <boost/type_traits.hpp>
#include <gtest/gtest.h>
#include <stan/math/prim/scal.hpp>

TEST(MetaTraits, containsFvar) {
  using stan::contains_fvar;
  EXPECT_FALSE(contains_fvar<double>::value);
}
