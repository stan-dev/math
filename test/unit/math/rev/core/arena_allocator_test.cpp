#include <stan/math/rev/core.hpp>
#include <gtest/gtest.h>
#include <vector>

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

  std::unordered_map<int, int, std::hash<int>, std::equal_to<int>,
                     stan::math::arena_allocator<std::pair<const int, int>>>
      m;
  m.insert(std::make_pair(3, 5));
  m.insert(std::make_pair(4, 6));
  for (auto it = m.begin(); it != m.end(); ++it) {
    EXPECT_TRUE(
        stan::math::ChainableStack::instance_->memalloc_.in_stack(&(*it)));
  }

  std::set<int, std::less<int>, stan::math::arena_allocator<int>> s;
  s.insert(3);
  s.insert(5);
  for (auto it = m.begin(); it != m.end(); ++it) {
    EXPECT_TRUE(
        stan::math::ChainableStack::instance_->memalloc_.in_stack(&(*it)));
  }
}

TEST(AgradRev, arena_allocator_test) {
  EXPECT_NO_THROW(arena_allocator_test());
}
