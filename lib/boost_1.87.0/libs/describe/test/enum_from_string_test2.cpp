// Copyright 2021 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/describe/enum_from_string.hpp>
#include <boost/describe/enum.hpp>
#include <boost/core/detail/string_view.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/config.hpp>

#if !defined(BOOST_DESCRIBE_CXX14)

#include <boost/config/pragma_message.hpp>

BOOST_PRAGMA_MESSAGE("Skipping test because C++14 is not available")
int main() {}

#else

enum E1 { v101 = 101, v102 = 102 };
BOOST_DESCRIBE_ENUM(E1, v101, v102)

enum class E2 { v201 = 201, v202 = 202 };
BOOST_DESCRIBE_ENUM(E2, v201, v202)

BOOST_DEFINE_ENUM(E3, v301, v302)
BOOST_DEFINE_ENUM_CLASS(E4, v401, v402)

template<class St> void test()
{
    using boost::describe::enum_from_string;

    {
        E1 w{};
        BOOST_TEST( enum_from_string( St( "v101" ), w ) ) && BOOST_TEST_EQ( w, v101 );
        BOOST_TEST( enum_from_string( St( "v102" ), w ) ) && BOOST_TEST_EQ( w, v102 );
        BOOST_TEST_NOT( enum_from_string( St( "v103" ), w ) );
    }

    {
        E2 w{};
        BOOST_TEST( enum_from_string( St( "v201" ), w ) ) && BOOST_TEST_EQ( (int)w, (int)E2::v201 );
        BOOST_TEST( enum_from_string( St( "v202" ), w ) ) && BOOST_TEST_EQ( (int)w, (int)E2::v202 );
        BOOST_TEST_NOT( enum_from_string( St( "v203" ), w ) );
    }

    {
        E3 w{};
        BOOST_TEST( enum_from_string( St( "v301" ), w ) ) && BOOST_TEST_EQ( w, v301 );
        BOOST_TEST( enum_from_string( St( "v302" ), w ) ) && BOOST_TEST_EQ( w, v302 );
        BOOST_TEST_NOT( enum_from_string( St( "v303" ), w ) );
    }

    {
        E4 w{};
        BOOST_TEST( enum_from_string( St( "v401" ), w ) ) && BOOST_TEST_EQ( (int)w, (int)E4::v401 );
        BOOST_TEST( enum_from_string( St( "v402" ), w ) ) && BOOST_TEST_EQ( (int)w, (int)E4::v402 );
        BOOST_TEST_NOT( enum_from_string( St( "v403" ), w ) );
    }
}

#if !defined(BOOST_NO_CXX17_HDR_STRING_VIEW)
# include <string_view>
#endif

int main()
{
    test<std::string>();
    test<boost::core::string_view>();

#if !defined(BOOST_NO_CXX17_HDR_STRING_VIEW)
    test<std::string_view>();
#endif

    return boost::report_errors();
}

#endif // !defined(BOOST_DESCRIBE_CXX14)
