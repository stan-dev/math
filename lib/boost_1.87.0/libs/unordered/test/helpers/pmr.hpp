
// Copyright 2024 Braden Ganetsky.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or move at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNORDERED_TEST_PMR_HEADER
#define BOOST_UNORDERED_TEST_PMR_HEADER

#include <boost/config.hpp>

#ifndef BOOST_NO_CXX17_HDR_MEMORY_RESOURCE

#include <memory_resource>

namespace test {
  class counted_new_delete_resource : public std::pmr::memory_resource
  {
  public:
    using std::pmr::memory_resource::memory_resource;
    ~counted_new_delete_resource() override {}

    std::size_t count() const { return _count; }

  private:
    void* do_allocate(std::size_t bytes, std::size_t alignment) override
    {
      _count += bytes;
      return ::operator new(bytes, std::align_val_t(alignment));
    }

    void do_deallocate(
      void* p, std::size_t bytes, std::size_t alignment) override
    {
      _count -= bytes;
#if __cpp_sized_deallocation
      ::operator delete(p, bytes, std::align_val_t(alignment));
#else
      ::operator delete(p, std::align_val_t(alignment));
#endif
    }

    bool do_is_equal(
      const std::pmr::memory_resource& other) const noexcept override
    {
      return this == &other;
    }

  private:
    std::size_t _count = 0;
  };
} // namespace test

#endif // !defined(BOOST_NO_CXX17_HDR_MEMORY_RESOURCE)

#endif // !defined(BOOST_UNORDERED_TEST_PMR_HEADER)
