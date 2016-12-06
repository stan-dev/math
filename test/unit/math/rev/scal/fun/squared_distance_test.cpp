#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/util.hpp>

TEST(MathRev, squared_distance) {
  double x1 = 1;
  double x2 = 4;
  stan::math::var v1, v2, f;
  std::vector<stan::math::var> vars;
  std::vector<double> grad_f;


  v1 = 1;
  v2 = 4;
  vars.push_back(v1);
  vars.push_back(v2);
  f = stan::math::squared_distance(v1, v2);
  f.grad(vars, grad_f);
  
  EXPECT_FLOAT_EQ(9, f.val());
  ASSERT_EQ(2, grad_f.size());
  EXPECT_FLOAT_EQ(-6, grad_f[0]);
  EXPECT_FLOAT_EQ(6, grad_f[1]);
  stan::math::recover_memory();
  vars.clear();


  v1 = 1;
  vars.push_back(v1);
  f = stan::math::squared_distance(v1, x2);
  f.grad(vars, grad_f);
  
  EXPECT_FLOAT_EQ(9, f.val());
  ASSERT_EQ(1, grad_f.size());
  EXPECT_FLOAT_EQ(-6, grad_f[0]);
  stan::math::recover_memory();
  vars.clear();


  
  v2 = 4;
  vars.push_back(v2);
  f = stan::math::squared_distance(x1, v2);
  f.grad(vars, grad_f);
  
  EXPECT_FLOAT_EQ(9, f.val());
  ASSERT_EQ(1, grad_f.size());
  EXPECT_FLOAT_EQ(6, grad_f[0]);
  stan::math::recover_memory();
  vars.clear();
}

TEST(MathRev, squared_distance_nan) {
  double x = 1;
  stan::math::var x_v = 1;
  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::math::var nan_v = std::numeric_limits<double>::quiet_NaN();
    
  EXPECT_THROW(stan::math::squared_distance(x_v, nan_v),
               std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(nan_v, x_v),
               std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(nan_v, nan_v),
               std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(x, nan_v),
               std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(nan, x_v),
               std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(nan, nan_v),
               std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(x_v, nan),
               std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(nan_v, x),
               std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(nan_v, nan),
               std::domain_error);
}


TEST(MathRev, squared_distance_inf) {
  double x = 1;
  stan::math::var x_v = 1;
  double inf = std::numeric_limits<double>::infinity();
  stan::math::var inf_v = std::numeric_limits<double>::infinity();
    
  EXPECT_THROW(stan::math::squared_distance(x_v, inf_v),
               std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(inf_v, x_v),
               std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(inf_v, inf_v),
               std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(x_v, inf),
               std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(inf_v, x),
               std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(inf_v, inf),
               std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(x, inf_v),
               std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(inf, x_v),
               std::domain_error);
  EXPECT_THROW(stan::math::squared_distance(inf, inf_v),
               std::domain_error);
}

TEST(MathRev, check_varis_on_stack) {
  stan::math::var v1 = 1;
  stan::math::var v2 = 4;
  test::check_varis_on_stack(stan::math::squared_distance(v1, v2));
}
