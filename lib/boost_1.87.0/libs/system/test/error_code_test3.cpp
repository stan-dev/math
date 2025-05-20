// Copyright 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt

#include <boost/system/error_code.hpp>
#include <boost/system/error_condition.hpp>
#include <boost/system/generic_category.hpp>
#include <boost/core/lightweight_test.hpp>

namespace sys = boost::system;

enum my_errc
{
    enomem_c = ENOMEM
};

enum my_errn
{
    enomem_n = ENOMEM
};

namespace boost {
namespace system {

template<> struct is_error_code_enum<my_errc>
{
    static const bool value = true;
};

template<> struct is_error_condition_enum<my_errn>
{
    static const bool value = true;
};

} // namespace system
} // namespace boost

sys::error_code make_error_code( my_errc e )
{
    return sys::error_code( e, sys::generic_category() );
}

sys::error_condition make_error_condition( my_errn e )
{
    return sys::error_condition( e, sys::generic_category() );
}

int main()
{
    sys::error_code ec = make_error_code( sys::errc::not_enough_memory );
    sys::error_condition en( sys::errc::not_enough_memory );

    BOOST_TEST_EQ( ec, en );
    BOOST_TEST_EQ( en, ec );

    BOOST_TEST_EQ( ec, enomem_c );
    BOOST_TEST_EQ( enomem_c, ec );

    BOOST_TEST_EQ( ec, enomem_n );
    BOOST_TEST_EQ( enomem_n, ec );

    BOOST_TEST_EQ( en, enomem_c );
    BOOST_TEST_EQ( enomem_c, en );

    BOOST_TEST_EQ( en, enomem_n );
    BOOST_TEST_EQ( enomem_n, en );

    return boost::report_errors();
}
