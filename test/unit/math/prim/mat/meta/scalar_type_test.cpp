#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <boost/type_traits.hpp>

TEST(MetaTraits, scalar_type) {
  using boost::is_same;
  using stan::scalar_type;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  stan::scalar_type<Matrix<double, Dynamic, Dynamic> >::type x = 1;
  EXPECT_FLOAT_EQ(1, x);

  // hack to get value of template into Google test macro 
  bool b3 = is_same<double, scalar_type<Matrix<double, Dynamic, Dynamic> >::type>::value;
  EXPECT_TRUE(b3);
}
