
// Copyright 2021-2023 Christian Mazakas.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/config.hpp>
#include <boost/config/pragma_message.hpp>
#include <boost/config/workaround.hpp>

#if BOOST_CXX_VERSION <= 199711L ||                                            \
  BOOST_WORKAROUND(BOOST_GCC_VERSION, < 40800) ||                              \
  BOOST_WORKAROUND(BOOST_MSVC, == 1900)

BOOST_PRAGMA_MESSAGE(
  "scoped allocator adaptor tests only work under C++11 and above")
int main() {}

#else

// This test is based on a user-submitted issue found here:
// https://github.com/boostorg/unordered/issues/22
//

#include "../helpers/unordered.hpp"

#include "../helpers/test.hpp"

#include <boost/cstdint.hpp>

#include <boost/core/ignore_unused.hpp>

#include <scoped_allocator>
#include <string>
#include <utility>
#include <vector>

namespace test {
  template <class T> struct allocator
  {
    typedef T value_type;

    allocator() = delete;
    allocator(int) {}
    allocator(allocator const&) = default;
    allocator(allocator&&) = default;

    template <class U> allocator(allocator<U> const&) {}

    BOOST_ATTRIBUTE_NODISCARD T* allocate(std::size_t n)
    {
      return static_cast<T*>(::operator new(n * sizeof(T)));
    }

    void deallocate(T* p, std::size_t) noexcept { ::operator delete(p); }

    bool operator==(allocator const&) const { return true; }
    bool operator!=(allocator const&) const { return false; }
  };
} // namespace test

typedef std::vector<boost::uint64_t, test::allocator<boost::uint64_t> >
  vector_type;

typedef std::pair<boost::uint64_t, vector_type> pair_type;

typedef std::scoped_allocator_adaptor<test::allocator<pair_type>,
  test::allocator<boost::uint64_t> >
  allocator_type;

template <class X> static void scoped_allocator(X*)
{
  allocator_type alloc(
    test::allocator<pair_type>(1337), test::allocator<boost::uint64_t>(7331));

  X map(alloc);

  for (unsigned i = 0; i < 10; ++i) {
    boost::ignore_unused(map[i]);
  }

  BOOST_TEST(map.size() == 10);
}

#ifdef BOOST_UNORDERED_FOA_TESTS
static boost::unordered_flat_map<const boost::uint64_t, vector_type,
  boost::hash<boost::uint64_t>, std::equal_to<boost::uint64_t>, allocator_type>*
  test_map;

static boost::unordered_node_map<const boost::uint64_t, vector_type,
  boost::hash<boost::uint64_t>, std::equal_to<boost::uint64_t>, allocator_type>*
  test_node_map;

UNORDERED_TEST(scoped_allocator, ((test_map)(test_node_map)))
#else
static boost::unordered_map<const boost::uint64_t, vector_type,
  boost::hash<boost::uint64_t>, std::equal_to<boost::uint64_t>, allocator_type>*
  test_map;

UNORDERED_TEST(scoped_allocator, ((test_map)))
#endif

RUN_TESTS()

#endif
