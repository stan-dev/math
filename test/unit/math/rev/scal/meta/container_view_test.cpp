#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <boost/type_traits/conditional.hpp>

TEST(MathMeta, container_view) {
  using stan::math::container_view;
  using stan::math::var;

  stan::math::start_nested();
  double y[1];
  container_view<var, double> view_test(var(4.0), y);
  view_test[0] = 1.0;
  EXPECT_FLOAT_EQ(1.0, view_test[0]);
  EXPECT_FLOAT_EQ(1.0, y[0]);

  var* test
    = static_cast<var*>(stan::math::vari::operator new(sizeof(var) * 1));
  container_view<var, var> view_test_v(var(4.0), test);
  view_test_v[0] = 2.0;
  EXPECT_FLOAT_EQ(2.0,view_test_v[0].val());
  EXPECT_FLOAT_EQ(2.0,test[0].val());
}

TEST(MathMeta, container_view_no_throw) {
  using stan::math::container_view;
  using boost::conditional;
  using stan::math::dummy;
  using stan::is_constant_struct;
  using stan::math::var;

  double arr[1];
  container_view<conditional<is_constant_struct<var>::value,dummy,var>::type, double> view_test(var(4.0), arr);
  EXPECT_NO_THROW(view_test[0]);
}
