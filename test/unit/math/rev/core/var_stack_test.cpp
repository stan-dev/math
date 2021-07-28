#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <gtest/gtest.h>
#include <vector>

struct foo : public stan::math::chainable_alloc {
  std::vector<double> x_;
  foo() : x_(3) {}
  ~foo() {}
  void chain() {}
};

TEST(AgradRev, varStackRecoverNestedSegFaultFix) {
  // this test failed in 2.5, but passes in 2.6
  stan::math::start_nested();
  foo* x = new foo();
  x->chain();
  stan::math::recover_memory_nested();
  // should be able to do this redundantly:
  stan::math::recover_memory();
}

// just test basic autodiff;  no more free_memory operation
TEST(AgradRev, varStack) {
  stan::math::var a = 2.0;
  stan::math::var b = -3.0;
  stan::math::var f = a * b;
  EXPECT_FLOAT_EQ(-6.0, f.val());

  std::vector<stan::math::var> x{a, b};
  std::vector<double> grad_f;
  f.grad(x, grad_f);
  EXPECT_FLOAT_EQ(-3.0, grad_f[0]);
  EXPECT_FLOAT_EQ(2.0, grad_f[1]);

  stan::math::var aa = 2.0;
  stan::math::var bb = -3.0;
  stan::math::var ff = aa * bb;
  EXPECT_FLOAT_EQ(-6.0, ff.val());

  std::vector<stan::math::var> xx = {aa, bb};
  std::vector<double> grad_ff;
  ff.grad(xx, grad_ff);
  EXPECT_FLOAT_EQ(-3.0, grad_ff[0]);
  EXPECT_FLOAT_EQ(2.0, grad_ff[1]);
}

TEST(AgradRev, recoverMemoryLogicError) {
  stan::math::start_nested();
  EXPECT_THROW(stan::math::recover_memory(), std::logic_error);
  // clean up for next test
  stan::math::recover_memory_nested();
}

TEST(AgradRev, recoverMemoryNestedLogicError) {
  EXPECT_THROW(stan::math::recover_memory_nested(), std::logic_error);
}
