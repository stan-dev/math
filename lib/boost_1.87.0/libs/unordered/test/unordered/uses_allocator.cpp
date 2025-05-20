
// Copyright 2023 Christian Mazakas.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "../helpers/test.hpp"

#include <boost/unordered/detail/implementation.hpp>
#include <boost/unordered/unordered_flat_map.hpp>
#include <boost/unordered/unordered_flat_set.hpp>
#include <boost/unordered/unordered_node_map.hpp>
#include <boost/unordered/unordered_node_set.hpp>

#include <boost/config.hpp>
#include <boost/config/pragma_message.hpp>
#include <boost/config/workaround.hpp>

#if BOOST_CXX_VERSION <= 199711L ||                                            \
  BOOST_WORKAROUND(BOOST_GCC_VERSION, < 40800) ||                              \
  (defined(BOOST_LIBSTDCXX_VERSION) && BOOST_CXX_VERSION > 201703L) ||         \
  (defined(BOOST_MSVC_FULL_VER) && BOOST_MSVC_FULL_VER >= 192000000 &&         \
    BOOST_MSVC_FULL_VER < 193000000)

// automatically disable this test for C++03 builds so we can use the STL's
// scoped_allocator_adaptor
// we remove C++20 support for libstdc++ builds because of:
// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=108952
//
// msvc-14.2 w/ C++20 is similarly affected
//


BOOST_PRAGMA_MESSAGE("uses_allocator tests require C++11, scoped_allocator")

int main() {}

#else

#include <memory>
#include <scoped_allocator>
#include <unordered_map>
#include <vector>

template <class T> struct allocator
{
  typedef T value_type;

  int tag_ = -1;

  allocator() = default;
  allocator(int tag) : tag_{tag} {}
  allocator(allocator const&) = default;
  allocator(allocator&&) = default;

  template <class U> allocator(allocator<U> const& rhs) : tag_{rhs.tag_} {}

  BOOST_ATTRIBUTE_NODISCARD T* allocate(std::size_t n)
  {
    return static_cast<T*>(::operator new(n * sizeof(T)));
  }

  void deallocate(T* p, std::size_t) noexcept { ::operator delete(p); }

  allocator& operator=(allocator const& rhs)
  {
    tag_ = rhs.tag_;
    return *this;
  }

  allocator& operator=(allocator&& rhs) noexcept
  {
    tag_ = rhs.tag_;
    return *this;
  }

  bool operator==(allocator const&) const { return true; }
  bool operator!=(allocator const&) const { return false; }
};

struct raii_tracker
{
  static int count;
  static int copy_count;
  static int move_count;
  static int alloc_move_count;

  using allocator_type = allocator<int>;

  allocator_type a_;

  raii_tracker(allocator_type a) : a_(a) { ++count; }
  raii_tracker(int, allocator_type const& a) : a_(a) { ++count; }

  raii_tracker(raii_tracker const&) { ++copy_count; }
  raii_tracker(raii_tracker&&) noexcept { ++move_count; }
  raii_tracker(raii_tracker&&, allocator_type const& a) noexcept : a_(a)
  {
    ++alloc_move_count;
  }

  allocator_type get_allocator() const noexcept { return a_; }

  friend bool operator==(raii_tracker const&, raii_tracker const&)
  {
    return true;
  }
};

int raii_tracker::count = 0;
int raii_tracker::copy_count = 0;
int raii_tracker::move_count = 0;
int raii_tracker::alloc_move_count = 0;

static void reset_counts()
{
  raii_tracker::count = 0;
  raii_tracker::copy_count = 0;
  raii_tracker::move_count = 0;
  raii_tracker::alloc_move_count = 0;
}

std::size_t hash_value(raii_tracker const&) { return 0; }

using map_allocator_type = std::scoped_allocator_adaptor<
  allocator<std::pair<raii_tracker const, raii_tracker> >, allocator<int> >;

using set_allocator_type =
  std::scoped_allocator_adaptor<allocator<raii_tracker>, allocator<int> >;

using map_type = boost::unordered_flat_map<raii_tracker, raii_tracker,
  boost::hash<raii_tracker>, std::equal_to<raii_tracker>, map_allocator_type>;

using node_map_type = boost::unordered_node_map<raii_tracker, raii_tracker,
  boost::hash<raii_tracker>, std::equal_to<raii_tracker>, map_allocator_type>;

using set_type = boost::unordered_flat_set<raii_tracker,
  boost::hash<raii_tracker>, std::equal_to<raii_tracker>, set_allocator_type>;

using node_set_type = boost::unordered_node_set<raii_tracker,
  boost::hash<raii_tracker>, std::equal_to<raii_tracker>, set_allocator_type>;

map_type* flat_map;
node_map_type* node_map;

set_type* flat_set;
node_set_type* node_set;

template <class X> static void map_uses_allocator_construction(X*)
{
  reset_counts();

  map_allocator_type alloc(
    allocator<std::pair<raii_tracker const, raii_tracker> >{12},
    allocator<int>{34});

  X map(1, alloc);
  map.emplace(
    std::piecewise_construct, std::make_tuple(1337), std::make_tuple(7331));

  BOOST_TEST_EQ(raii_tracker::count, 2);
  BOOST_TEST_EQ(raii_tracker::move_count, 0);
  BOOST_TEST_EQ(raii_tracker::alloc_move_count, 2);

  BOOST_TEST_EQ(map.begin()->first.get_allocator().tag_, 34);
  BOOST_TEST_EQ(map.begin()->second.get_allocator().tag_, 34);
}

template <class X> static void set_uses_allocator_construction(X*)
{
  reset_counts();

  set_allocator_type alloc(allocator<raii_tracker>{12}, allocator<int>{34});

  X set(1, alloc);
  set.emplace();

  BOOST_TEST_EQ(raii_tracker::count, 1);
  BOOST_TEST_EQ(raii_tracker::move_count, 0);
  BOOST_TEST_EQ(raii_tracker::alloc_move_count, 1);

  BOOST_TEST_EQ(set.begin()->get_allocator().tag_, 34);
}

UNORDERED_TEST(map_uses_allocator_construction, ((flat_map)(node_map)))
UNORDERED_TEST(set_uses_allocator_construction, ((flat_set)(node_set)))

RUN_TESTS()

#endif
