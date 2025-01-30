//////////////////////////////////////////////////////////////////////////////
//
// (C) Copyright David Abrahams, Vicente Botet, Ion Gaztanaga 2009.
// (C) Copyright Ion Gaztanaga 2009-2014.
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/move for documentation.
//
//////////////////////////////////////////////////////////////////////////////
#include "../example/movable.hpp"
#include "../example/copymovable.hpp"

#ifdef BOOST_MOVE_MOVE_UTILITY_CORE_HPP
#error "<boost/move/utility_core.hpp> should not be included for this test. Minimal headers are required."
#endif

movable produce_movable()
{  return movable(); }

movable produce_copyable()
{  return movable(); }

int main()
{
   movable m;
   copyable c;
   ::boost::movelib::ignore(m);
   ::boost::movelib::ignore(c);

   movable m2 = produce_movable();
   movable c2 = produce_copyable();

   ::boost::movelib::ignore(m2);
   ::boost::movelib::ignore(c2);

   return 0;
}
