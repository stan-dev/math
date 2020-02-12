#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/core/gradable.hpp>
#include <gtest/gtest.h>
#include <vector>

struct AgradLocalNested : public testing::Test {
  void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
  }
};

TEST_F(AgradLocalNested, local_nested_autodiff_base) {
  {
    stan::math::local_nested_autodiff nested;
    EXPECT_THROW(stan::math::recover_memory(), std::logic_error);
  }
  stan::math::recover_memory();  // Should not throw

  gradable g_out = setup_quad_form();
  for (int i = 0; i < 100; ++i) {
    stan::math::local_nested_autodiff nested;
    gradable g = setup_quad_form();
    g.test();
    nested.set_zero_all_adjoints();
    EXPECT_EQ(g.adj(), 0);
  }
  g_out.test();
}

TEST_F(AgradLocalNested, local_nested_autodiff_Gradient1) {
  using stan::math::local_nested_autodiff;

  gradable g0 = setup_simple();

  {
    local_nested_autodiff nested;
    gradable g1 = setup_quad_form();
    g1.test();
  }

  {
    local_nested_autodiff nested;
    gradable g2 = setup_simple();
    g2.test();
  }

  g0.test();
  stan::math::recover_memory();
}

TEST_F(AgradLocalNested, local_nested_autodiff_Gradient2) {
  using stan::math::local_nested_autodiff;

  gradable g0 = setup_quad_form();

  {
    local_nested_autodiff nested;
    gradable g1 = setup_simple();
    g1.test();
  }

  {
    local_nested_autodiff nested;
    gradable g2 = setup_quad_form();
    g2.test();
  }

  g0.test();
  stan::math::recover_memory();
}

TEST_F(AgradLocalNested, local_nested_autodiff_Gradient3) {
  using stan::math::local_nested_autodiff;

  {
    local_nested_autodiff nested;
    gradable g1 = setup_simple();
    {
      local_nested_autodiff nested2;
      gradable g2 = setup_quad_form();
      {
        local_nested_autodiff nested3;
        gradable g3 = setup_quad_form();
        {
          local_nested_autodiff nested4;
          gradable g4 = setup_simple();
          g4.test();
        }
        g3.test();
      }
      g2.test();
    }
    g1.test();
  }
  stan::math::recover_memory();
}