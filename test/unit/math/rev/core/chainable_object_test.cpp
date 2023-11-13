#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

class ChainableObjectTest {
 public:
  static int counter;

  ~ChainableObjectTest() { counter++; }
};

int ChainableObjectTest::counter = 0;

TEST(AgradRevChain, chainable_object_test) {
  {
    EXPECT_NO_THROW(new stan::math::chainable_object<ChainableObjectTest>(
        ChainableObjectTest()));
    ChainableObjectTest::counter = 0;
  }

  EXPECT_EQ((ChainableObjectTest::counter), 0);
  stan::math::recover_memory();
  EXPECT_EQ((ChainableObjectTest::counter), 1);
}

TEST(AgradRevChain, chainable_object_nested_test) {
  stan::math::start_nested();

  {
    EXPECT_NO_THROW(new stan::math::chainable_object<ChainableObjectTest>(
        ChainableObjectTest()));
    ChainableObjectTest::counter = 0;
  }

  EXPECT_EQ((ChainableObjectTest::counter), 0);

  stan::math::recover_memory_nested();

  EXPECT_EQ((ChainableObjectTest::counter), 1);
}

TEST(AgradRevChain, make_chainable_ptr_test) {
  {
    EXPECT_NO_THROW(stan::math::make_chainable_ptr(ChainableObjectTest()));
    ChainableObjectTest::counter = 0;
  }

  EXPECT_EQ((ChainableObjectTest::counter), 0);
  stan::math::recover_memory();
  EXPECT_EQ((ChainableObjectTest::counter), 1);
}

TEST(AgradRevChain, make_chainable_ptr_nested_test) {
  stan::math::start_nested();

  {
    EXPECT_NO_THROW(stan::math::make_chainable_ptr(ChainableObjectTest()));
    ChainableObjectTest::counter = 0;
  }

  EXPECT_EQ((ChainableObjectTest::counter), 0);

  stan::math::recover_memory_nested();

  EXPECT_EQ((ChainableObjectTest::counter), 1);
}

class UnsafeChainableObjectTest {
 public:
  static int counter;

  ~UnsafeChainableObjectTest() { counter++; }
};

int UnsafeChainableObjectTest::counter = 0;

TEST(AgradRevChain, unsafe_chainable_object_test) {
  {
    EXPECT_NO_THROW(
        new stan::math::unsafe_chainable_object<UnsafeChainableObjectTest>(
            UnsafeChainableObjectTest()));
    UnsafeChainableObjectTest::counter = 0;
  }

  EXPECT_EQ((UnsafeChainableObjectTest::counter), 0);
  stan::math::recover_memory();
  EXPECT_EQ((UnsafeChainableObjectTest::counter), 1);
}

TEST(AgradRevChain, unsafe_chainable_object_nested_test) {
  stan::math::start_nested();

  {
    EXPECT_NO_THROW(
        stan::math::make_unsafe_chainable_ptr(UnsafeChainableObjectTest()));
    UnsafeChainableObjectTest::counter = 0;
  }

  EXPECT_EQ((UnsafeChainableObjectTest::counter), 0);

  stan::math::recover_memory_nested();

  EXPECT_EQ((UnsafeChainableObjectTest::counter), 1);
}

TEST(AgradRevChain, make_unsafe_chainable_ptr_test) {
  {
    EXPECT_NO_THROW(
        stan::math::make_unsafe_chainable_ptr(UnsafeChainableObjectTest()));
    UnsafeChainableObjectTest::counter = 0;
  }

  EXPECT_EQ((UnsafeChainableObjectTest::counter), 0);
  stan::math::recover_memory();
  EXPECT_EQ((UnsafeChainableObjectTest::counter), 1);
}

TEST(AgradRevChain, make_unsafe_chainable_ptr_nested_test) {
  stan::math::start_nested();

  {
    EXPECT_NO_THROW(
        stan::math::make_unsafe_chainable_ptr(UnsafeChainableObjectTest()));
    UnsafeChainableObjectTest::counter = 0;
  }

  EXPECT_EQ((UnsafeChainableObjectTest::counter), 0);

  stan::math::recover_memory_nested();

  EXPECT_EQ((UnsafeChainableObjectTest::counter), 1);
}
