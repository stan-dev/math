#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

class ChainableObjectTest {
 public:
  static int counter;

  ~ChainableObjectTest() { counter++; }
};

int ChainableObjectTest::counter = 0;

TEST(AgradRev, chainable_object_test) {
  {
    auto ptr = new stan::math::chainable_object<ChainableObjectTest>(
        ChainableObjectTest());
    ChainableObjectTest::counter = 0;
  }

  EXPECT_EQ((ChainableObjectTest::counter), 0);
  stan::math::recover_memory();
  EXPECT_EQ((ChainableObjectTest::counter), 1);
}

TEST(AgradRev, chainable_object_nested_test) {
  stan::math::start_nested();

  {
    auto ptr = new stan::math::chainable_object<ChainableObjectTest>(
        ChainableObjectTest());
    ChainableObjectTest::counter = 0;
  }

  EXPECT_EQ((ChainableObjectTest::counter), 0);

  stan::math::recover_memory_nested();

  EXPECT_EQ((ChainableObjectTest::counter), 1);
}

TEST(AgradRev, make_chainable_ptr_test) {
  {
    ChainableObjectTest* ptr
        = stan::math::make_chainable_ptr(ChainableObjectTest());
    ChainableObjectTest::counter = 0;
  }

  EXPECT_EQ((ChainableObjectTest::counter), 0);
  stan::math::recover_memory();
  EXPECT_EQ((ChainableObjectTest::counter), 1);
}

TEST(AgradRev, make_chainable_ptr_nested_test) {
  stan::math::start_nested();

  {
    ChainableObjectTest* ptr
        = stan::math::make_chainable_ptr(ChainableObjectTest());
    ChainableObjectTest::counter = 0;
  }

  EXPECT_EQ((ChainableObjectTest::counter), 0);

  stan::math::recover_memory_nested();

  EXPECT_EQ((ChainableObjectTest::counter), 1);
}
