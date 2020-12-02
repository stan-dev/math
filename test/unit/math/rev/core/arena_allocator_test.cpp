#include <stan/math/rev/core.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <unordered_set>

void arena_allocator_test() {
  std::vector<int, stan::math::arena_allocator<int>> v;
  v.reserve(3);
  for (int i = 0; i < 5; i++) {
    v.push_back(1);
  }
  EXPECT_TRUE(
      stan::math::ChainableStack::instance_->memalloc_.in_stack(v.data()));
  EXPECT_TRUE(stan::math::ChainableStack::instance_->memalloc_.in_stack(
      v.data() + v.size()));

  std::vector<int, stan::math::arena_allocator<int>> v2(2, 3);
  for (int i = 0; i < 5; i++) {
    v2.push_back(2);
  }
  EXPECT_TRUE(
      stan::math::ChainableStack::instance_->memalloc_.in_stack(v2.data()));
  EXPECT_TRUE(stan::math::ChainableStack::instance_->memalloc_.in_stack(
      v2.data() + v2.size()));

  std::vector<int, stan::math::arena_allocator<int>> v3(v2);
  for (int i = 0; i < 15; i++) {
    v3.push_back(3);
  }
  EXPECT_TRUE(
      stan::math::ChainableStack::instance_->memalloc_.in_stack(v3.data()));
  EXPECT_TRUE(stan::math::ChainableStack::instance_->memalloc_.in_stack(
      v3.data() + v3.size()));
}

TEST(AgradRev, arena_allocator_test) {
  EXPECT_NO_THROW(arena_allocator_test());
}

TEST(AgradRev, arena_allocator_unorderedset_test) {
  std::unordered_set<int, std::hash<int>, std::equal_to<int>,
                     stan::math::arena_allocator<int>>
      x_test;
  x_test.reserve(10);
  for (int i = 0; i < 5; i++) {
    x_test.insert(i);
    x_test.insert(i);
  }
  auto x_iter = x_test.begin();
  while (x_iter != x_test.end()) {
    EXPECT_TRUE(
        stan::math::ChainableStack::instance_->memalloc_.in_stack(&(*x_iter)));
    x_iter++;
  }
  auto x_test2 = x_test;
  auto x_iter2 = x_test2.begin();
  while (x_iter2 != x_test2.end()) {
    EXPECT_TRUE(
        stan::math::ChainableStack::instance_->memalloc_.in_stack(&(*x_iter2)));
    x_iter2++;
  }
}
