// Copyright 2024 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/smart_ptr/detail/sp_type_traits.hpp>
#include <boost/core/lightweight_test.hpp>

int main()
{
    using boost::detail::sp_type_with_alignment;

    BOOST_TEST_EQ( alignof( sp_type_with_alignment<1>::type ), 1 );
    BOOST_TEST_EQ( alignof( sp_type_with_alignment<2>::type ), 2 );
    BOOST_TEST_EQ( alignof( sp_type_with_alignment<4>::type ), 4 );
    BOOST_TEST_EQ( alignof( sp_type_with_alignment<8>::type ), 8 );

    return boost::report_errors();
}
