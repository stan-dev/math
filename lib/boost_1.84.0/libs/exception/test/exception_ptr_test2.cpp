// Copyright 2022 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/config.hpp>

#if defined(BOOST_NO_EXCEPTIONS)

#include <boost/config/pragma_message.hpp>
BOOST_PRAGMA_MESSAGE( "Skipping test because BOOST_NO_EXCEPTIONS is defined" )

int main() {}

#elif defined(BOOST_NO_CXX11_HDR_EXCEPTION)

#include <boost/config/pragma_message.hpp>
BOOST_PRAGMA_MESSAGE( "Skipping test because BOOST_NO_CXX11_HDR_EXCEPTION is defined" )

int main() {}

#else

#include <boost/exception_ptr.hpp>
#include <boost/exception/exception.hpp>
#include <boost/core/lightweight_test.hpp>
#include <exception>
#include <new>
#include <stdexcept>

class my_exception
{
};

class my_exception2: public std::exception
{
};

class my_exception3: public std::bad_alloc
{
};

class my_exception4: public std::exception, public boost::exception
{
};

class my_exception5: public std::logic_error, public virtual boost::exception
{
public:

    my_exception5(): std::logic_error( "" ) {}
};

int main()
{
    try
    {
        throw my_exception();
    }
    catch( ... )
    {
        boost::exception_ptr p = boost::current_exception();
        BOOST_TEST_THROWS( boost::rethrow_exception( p ), my_exception );
    }

    try
    {
        throw my_exception2();
    }
    catch( ... )
    {
        boost::exception_ptr p = boost::current_exception();
        BOOST_TEST_THROWS( boost::rethrow_exception( p ), my_exception2 );
        BOOST_TEST_THROWS( boost::rethrow_exception( p ), std::exception );
    }

    try
    {
        throw my_exception3();
    }
    catch( ... )
    {
        boost::exception_ptr p = boost::current_exception();
        BOOST_TEST_THROWS( boost::rethrow_exception( p ), my_exception3 );
        BOOST_TEST_THROWS( boost::rethrow_exception( p ), std::bad_alloc );
    }

    try
    {
        throw my_exception4();
    }
    catch( ... )
    {
        boost::exception_ptr p = boost::current_exception();
        BOOST_TEST_THROWS( boost::rethrow_exception( p ), my_exception4 );
        BOOST_TEST_THROWS( boost::rethrow_exception( p ), std::exception );
        BOOST_TEST_THROWS( boost::rethrow_exception( p ), boost::exception );
    }

    try
    {
        throw my_exception5();
    }
    catch( ... )
    {
        boost::exception_ptr p = boost::current_exception();
        BOOST_TEST_THROWS( boost::rethrow_exception( p ), my_exception5 );
        BOOST_TEST_THROWS( boost::rethrow_exception( p ), std::logic_error );
        BOOST_TEST_THROWS( boost::rethrow_exception( p ), boost::exception );
    }

    return boost::report_errors();
}

#endif
