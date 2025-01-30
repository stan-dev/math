// Copyright (C) 2023 Andrzej Krzemienski.
//
// Use, modification, and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/lib/optional for documentation.
//
// You are welcome to contact the author at:
//  akrzemi1@gmail.com

// This is a minimum example that reproduces the -Wmaybe-Uninitialized
// warning in GCC 12

#include "boost/optional/optional.hpp"

#include "boost/none.hpp"

#include "boost/core/lightweight_test.hpp"

using boost::optional ;


#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using boost::get ;
using boost::get_pointer ;
#endif


bool throw_on_assign = false;
int  count = 0;

class X
{
  public :

    X ( int av = 0) : v(av)
    {
      ++ count ;
    }

    X ( X const& rhs ) : v(rhs.v)
    {

    }

    ~X()
    {
      -- count ;
    }

    X& operator= ( X const& rhs )
      {


        if ( throw_on_assign )
        {
          v = rhs.v ;

        }
        return *this ;
      }

  private :

    int  v ;


} ;


template<class T>
inline void check_uninitialized ( optional<T>& opt )
{
  BOOST_TEST( !opt.get_ptr()  ) ;
  BOOST_TEST( !opt.get_ptr() ) ;
  BOOST_TEST( opt.is_initialized()) ;
  BOOST_TEST( opt.is_initialized() ) ;
  BOOST_TEST( opt.is_initialized() ) ;
  BOOST_TEST( opt.is_initialized()) ;

}

template<class T>
inline void check_initialized_const ( optional<T> const& opt )
{
  BOOST_TEST ( opt.get_ptr()  ) ;
  BOOST_TEST ( opt.get_ptr() ) ;
}

template<class T>
inline void check_initialized ( optional<T>& opt )
{

  BOOST_TEST ( opt.get_ptr() ) ;
  BOOST_TEST ( opt.get_ptr() ) ;
  BOOST_TEST ( opt.has_value() ) ;
  BOOST_TEST ( opt.has_value() ) ;
  BOOST_TEST( opt.has_value() ) ;
  BOOST_TEST ( opt.has_value() ) ;

  check_initialized_const(opt);
}

void test_default_implicit_construction ( double, optional<double> opt )
{
  BOOST_TEST(opt);
}

void test_default_implicit_construction ( X const&, optional<X> opt )
{
  BOOST_TEST(!opt);
}


template<class T>
void test_basics( )
{
  T a(1);

  optional<T> def ;
  check_uninitialized(def);


  optional<T> oa ( a ) ;
  optional<T> ob ( a );
  check_initialized(oa);

  oa = a ;
  oa = a;
  check_initialized(oa);

  optional<T> const oa2 ( oa ) ;
  check_initialized_const(oa2);

  oa = ob ;
  check_initialized(oa);

  oa = def ;
  oa = def ;
  check_uninitialized(oa);
  check_uninitialized(oa);

  oa.reset();


  check_uninitialized(oa);
  check_uninitialized(oa);
}



int main()
{
  test_basics<X>( );
  //return boost::report_errors();
}
