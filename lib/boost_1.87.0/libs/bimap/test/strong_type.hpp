// Boost.Bimap
//
// Copyright (c) 2024 Joaquin M Lopez Munoz
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef LIBS_BIMAP_TEST_STRONG_TYPE_HPP
#define LIBS_BIMAP_TEST_STRONG_TYPE_HPP

#if defined(_MSC_VER)
#pragma once
#endif

#include <boost/container_hash/hash.hpp>

// std
#include <cstddef>
#include <stdexcept>

template< class T >
struct strong
{
  template< class Q >
  strong(const Q& x_):x(x_){}

  T x;
};

template< class T >
bool operator<(const strong<T>& x,const strong<T>& y)
{
  return x.x<y.x;
}

template< class T >
bool operator==(const strong<T>& x,const strong<T>& y)
{
  return x.x==y.x;
}

template< class T >
std::size_t hash_value(const strong<T>& x)
{
  return boost::hash<T>()(x.x);
}

template< class T >
struct semistrong: strong<T>
{
  using strong<T>::strong;

  // semistrong<T> is formally convertible to T but throws when
  // conversion actually called

  operator const T&() const
  {
    throw std::runtime_error("semistrong<T>  -> T conversion called");
  }
};

#endif // LIBS_BIMAP_TEST_STRONG_TYPE_HPP
