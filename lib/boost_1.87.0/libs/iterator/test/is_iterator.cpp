/*
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * https://www.boost.org/LICENSE_1_0.txt)
 *
 * Copyright (c) 2023 Andrey Semashev
 */
/*!
 * \file is_iterator.cpp
 *
 * This header contains tests for the \c is_iterator type trait.
 */

#include <boost/iterator/is_iterator.hpp>
#include <cstddef>
#include <list>
#include <vector>
#include <string>
#include <iterator>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/core/lightweight_test.hpp>

template< typename T >
struct value_iterator
{
    typedef std::input_iterator_tag iterator_category;
    typedef T value_type;
    typedef std::ptrdiff_t difference_type;
    typedef T* pointer;
    typedef T& reference;

    value_type operator*() const;
};

template< typename T >
struct proxy_iterator
{
    typedef T value_type;
    typedef std::output_iterator_tag iterator_category;
    typedef std::ptrdiff_t difference_type;
    typedef T* pointer;
    typedef T& reference;

    struct proxy
    {
        operator value_type&() const;
        proxy& operator=(value_type) const;
    };

    proxy operator*() const;
};

template< typename T >
struct lvalue_iterator
{
    typedef T value_type;
    typedef T& reference;
    typedef T difference_type;
    typedef std::input_iterator_tag iterator_category;
    typedef T* pointer;

    T& operator*() const;
    lvalue_iterator& operator++();
    lvalue_iterator operator++(int);
};

template< typename T >
struct constant_lvalue_iterator
{
    typedef T value_type;
    typedef T const& reference;
    typedef T difference_type;
    typedef std::input_iterator_tag iterator_category;
    typedef T const* pointer;

    T const& operator*() const;
    constant_lvalue_iterator& operator++();
    constant_lvalue_iterator operator++(int);
};

template< typename Iterator >
class adapted_iterator :
    public boost::iterators::iterator_adaptor< adapted_iterator< Iterator >, Iterator >
{
    friend class iterator_core_access;

private:
    typedef boost::iterators::iterator_adaptor< adapted_iterator< Iterator >, Iterator > base_type;

private:
    typename base_type::reference dereference() const;
    void increment();
    void decrement();
    void advance(typename base_type::difference_type n);
    template< typename OtherIterator >
    typename base_type::difference_type distance_to(adapted_iterator< OtherIterator > const& y) const;
};

struct complete {};
struct incomplete;

int main()
{
    BOOST_TEST(boost::iterators::is_iterator< int* >::value);
    BOOST_TEST(boost::iterators::is_iterator< const int* >::value);
    BOOST_TEST(boost::iterators::is_iterator< complete* >::value);
    BOOST_TEST(boost::iterators::is_iterator< std::reverse_iterator< int* > >::value);
    BOOST_TEST(boost::iterators::is_iterator< std::reverse_iterator< complete* > >::value);
    BOOST_TEST(boost::iterators::is_iterator< adapted_iterator< int* > >::value);

    BOOST_TEST(boost::iterators::is_iterator< std::string::iterator >::value);
    BOOST_TEST(boost::iterators::is_iterator< std::string::const_iterator >::value);
    BOOST_TEST(boost::iterators::is_iterator< std::string::reverse_iterator >::value);
    BOOST_TEST(boost::iterators::is_iterator< std::string::const_reverse_iterator >::value);

    BOOST_TEST(boost::iterators::is_iterator< std::list< int >::iterator >::value);
    BOOST_TEST(boost::iterators::is_iterator< std::list< int >::const_iterator >::value);
    BOOST_TEST(boost::iterators::is_iterator< std::list< int >::reverse_iterator >::value);
    BOOST_TEST(boost::iterators::is_iterator< std::list< int >::const_reverse_iterator >::value);

    BOOST_TEST(boost::iterators::is_iterator< std::vector< int >::iterator >::value);
    BOOST_TEST(boost::iterators::is_iterator< std::vector< int >::const_iterator >::value);
    BOOST_TEST(boost::iterators::is_iterator< std::vector< int >::reverse_iterator >::value);
    BOOST_TEST(boost::iterators::is_iterator< std::vector< int >::const_reverse_iterator >::value);

    BOOST_TEST(boost::iterators::is_iterator< std::insert_iterator< std::vector< int > > >::value);
    BOOST_TEST(boost::iterators::is_iterator< std::back_insert_iterator< std::vector< int > > >::value);
    BOOST_TEST(boost::iterators::is_iterator< std::front_insert_iterator< std::vector< int > > >::value);
    BOOST_TEST(boost::iterators::is_iterator< std::istream_iterator< int > >::value);
    BOOST_TEST(boost::iterators::is_iterator< std::ostream_iterator< int > >::value);
    BOOST_TEST(boost::iterators::is_iterator< std::istreambuf_iterator< char > >::value);
    BOOST_TEST(boost::iterators::is_iterator< std::ostreambuf_iterator< char > >::value);

    BOOST_TEST(!boost::iterators::is_iterator< int >::value);
    BOOST_TEST(!boost::iterators::is_iterator< complete >::value);
    BOOST_TEST(!boost::iterators::is_iterator< void >::value);
    BOOST_TEST(!boost::iterators::is_iterator< const void >::value);
    BOOST_TEST(!boost::iterators::is_iterator< void* >::value);
#if defined(BOOST_TT_HAS_WORKING_IS_COMPLETE)
    BOOST_TEST(!boost::iterators::is_iterator< incomplete >::value);
    BOOST_TEST(!boost::iterators::is_iterator< incomplete* >::value);
#endif
    BOOST_TEST(!boost::iterators::is_iterator< int (int) >::value);
    BOOST_TEST(!boost::iterators::is_iterator< int (*)(int) >::value);
    BOOST_TEST(!boost::iterators::is_iterator< int complete::* >::value);
    BOOST_TEST(!boost::iterators::is_iterator< int (complete::*)(int) >::value);
    BOOST_TEST(!boost::iterators::is_iterator< int (complete::*)(int) const >::value);
#if defined(__cpp_noexcept_function_type) && (__cpp_noexcept_function_type >= 201510L)
    BOOST_TEST(!boost::iterators::is_iterator< int (*)(int) noexcept >::value);
    BOOST_TEST(!boost::iterators::is_iterator< int (complete::*)(int) noexcept >::value);
    BOOST_TEST(!boost::iterators::is_iterator< int (complete::*)(int) const noexcept >::value);
#endif
    BOOST_TEST(!boost::iterators::is_iterator< int[] >::value);
    BOOST_TEST(!boost::iterators::is_iterator< int[10] >::value);
    BOOST_TEST(!boost::iterators::is_iterator< int*[] >::value);
    BOOST_TEST(!boost::iterators::is_iterator< int*[10] >::value);

    BOOST_TEST(!boost::iterators::is_iterator< int& >::value);
    BOOST_TEST(!boost::iterators::is_iterator< int*& >::value);
    BOOST_TEST(!boost::iterators::is_iterator< int (&)(int) >::value);
    BOOST_TEST(!boost::iterators::is_iterator< int (&)[10] >::value);

    return boost::report_errors();
}
