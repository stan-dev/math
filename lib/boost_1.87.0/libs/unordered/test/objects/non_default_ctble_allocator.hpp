
// Copyright 2024 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#if !defined(BOOST_UNORDERED_TEST_NON_DEFAULT_CTBLE_ALLOCATOR_HEADER)
#define BOOST_UNORDERED_TEST_NON_DEFAULT_CTBLE_ALLOCATOR_HEADER

#include <memory>

namespace test
{
  template <class T> struct non_default_ctble_allocator
  {
    typedef T value_type;

    non_default_ctble_allocator(int) {}

    template <class U> 
    non_default_ctble_allocator(const non_default_ctble_allocator<U>&) {}

    template<class U>
    bool operator==(non_default_ctble_allocator<U> const &) const noexcept
    {
      return true;
    }

    template<class U>
    bool operator!=(non_default_ctble_allocator<U> const &) const noexcept
    {
      return false;
    }
    
    T* allocate(std::size_t n)
    {
      return std::allocator<T>().allocate(n); 
    }

    void deallocate(T* p, std::size_t n)
    {
      std::allocator<T>().deallocate(p, n); 
    }
  };
}

#endif
