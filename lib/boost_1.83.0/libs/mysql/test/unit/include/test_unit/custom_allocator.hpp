//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CUSTOM_ALLOCATOR_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CUSTOM_ALLOCATOR_HPP

// This is used for concept checks. Not implemented anywhere.
#include <cstddef>

namespace boost {
namespace mysql {
namespace test {

template <class T>
struct custom_allocator
{
    using value_type = T;
    custom_allocator() noexcept;
    template <class U>
    custom_allocator(const custom_allocator<U>&) noexcept;
    T* allocate(std::size_t n);
    void deallocate(T* p, std::size_t n);
};

template <class T>
struct custom_allocator_no_defctor : custom_allocator<T>
{
    custom_allocator_no_defctor(int) noexcept {}
};

template <class T, class U>
constexpr bool operator==(const custom_allocator<T>&, const custom_allocator<U>&) noexcept;

template <class T, class U>
constexpr bool operator!=(const custom_allocator<T>&, const custom_allocator<U>&) noexcept;

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
