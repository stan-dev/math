#include <stan/math/rev/core.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <deque>
#include <list>
#include <forward_list>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>

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

  std::deque<int, stan::math::arena_allocator<int>> d(4, 1);
  for (auto it = d.begin(); it != d.end(); ++it) {
    EXPECT_TRUE(
        stan::math::ChainableStack::instance_->memalloc_.in_stack(&(*it)));
  }

  std::list<int, stan::math::arena_allocator<int>> l(4, 1);
  for (auto it = l.begin(); it != l.end(); ++it) {
    EXPECT_TRUE(
        stan::math::ChainableStack::instance_->memalloc_.in_stack(&(*it)));
  }

  std::forward_list<int, stan::math::arena_allocator<int>> fl(4, 1);
  for (auto it = fl.begin(); it != fl.end(); ++it) {
    EXPECT_TRUE(
        stan::math::ChainableStack::instance_->memalloc_.in_stack(&(*it)));
  }

  std::set<int, std::less<int>, stan::math::arena_allocator<int>> s;
  s.insert(3);
  s.insert(5);
  for (auto it = s.begin(); it != s.end(); ++it) {
    EXPECT_TRUE(
        stan::math::ChainableStack::instance_->memalloc_.in_stack(&(*it)));
  }

  std::multiset<int, std::less<int>, stan::math::arena_allocator<int>> ms;
  ms.insert(3);
  ms.insert(5);
  ms.insert(5);
  for (auto it = ms.begin(); it != ms.end(); ++it) {
    EXPECT_TRUE(
        stan::math::ChainableStack::instance_->memalloc_.in_stack(&(*it)));
  }

  std::unordered_set<int, std::hash<int>, std::equal_to<int>,
                     stan::math::arena_allocator<int>>
      u_s;
  u_s.insert(3);
  u_s.insert(5);
  for (auto it = u_s.begin(); it != u_s.end(); ++it) {
    EXPECT_TRUE(
        stan::math::ChainableStack::instance_->memalloc_.in_stack(&(*it)));
  }

  std::unordered_multiset<int, std::hash<int>, std::equal_to<int>,
                          stan::math::arena_allocator<int>>
      u_ms;
  u_ms.insert(3);
  u_ms.insert(5);
  u_ms.insert(5);
  for (auto it = u_ms.begin(); it != u_ms.end(); ++it) {
    EXPECT_TRUE(
        stan::math::ChainableStack::instance_->memalloc_.in_stack(&(*it)));
  }

  std::map<int, int, std::less<int>,
           stan::math::arena_allocator<std::pair<const int, int>>>
      m;
  m.insert(std::make_pair(3, 5));
  m.insert(std::make_pair(4, 6));
  for (auto it = m.begin(); it != m.end(); ++it) {
    EXPECT_TRUE(
        stan::math::ChainableStack::instance_->memalloc_.in_stack(&(*it)));
  }

  std::multimap<int, int, std::less<int>,
                stan::math::arena_allocator<std::pair<const int, int>>>
      mm;
  mm.insert(std::make_pair(3, 5));
  mm.insert(std::make_pair(4, 6));
  mm.insert(std::make_pair(4, 8));
  for (auto it = mm.begin(); it != mm.end(); ++it) {
    EXPECT_TRUE(
        stan::math::ChainableStack::instance_->memalloc_.in_stack(&(*it)));
  }

  std::unordered_map<int, int, std::hash<int>, std::equal_to<int>,
                     stan::math::arena_allocator<std::pair<const int, int>>>
      u_m;
  u_m.insert(std::make_pair(3, 5));
  u_m.insert(std::make_pair(4, 6));
  for (auto it = u_m.begin(); it != u_m.end(); ++it) {
    EXPECT_TRUE(
        stan::math::ChainableStack::instance_->memalloc_.in_stack(&(*it)));
  }

  std::unordered_multimap<
      int, int, std::hash<int>, std::equal_to<int>,
      stan::math::arena_allocator<std::pair<const int, int>>>
      u_mm;
  u_mm.insert(std::make_pair(3, 5));
  u_mm.insert(std::make_pair(4, 6));
  u_mm.insert(std::make_pair(4, 8));
  for (auto it = u_mm.begin(); it != u_mm.end(); ++it) {
    EXPECT_TRUE(
        stan::math::ChainableStack::instance_->memalloc_.in_stack(&(*it)));
  }
}

TEST(AgradRevArena, arena_allocator_test) {
  EXPECT_NO_THROW(arena_allocator_test());
}
