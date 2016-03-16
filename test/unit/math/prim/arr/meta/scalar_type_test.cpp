#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <boost/type_traits.hpp>

TEST(MetaTraits, scalar_type) {
  using boost::is_same;
  using stan::scalar_type;
  using std::vector;

  stan::scalar_type<std::vector<int> >::type n = 1;
  EXPECT_EQ(1,n);

  // hack to get value of template into Google test macro 
  bool b3 = is_same<double, scalar_type<vector<double> >::type>::value;
  EXPECT_TRUE(b3);
}
