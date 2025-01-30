// Copyright 2023 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt

#include <boost/system/error_category.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/core/snprintf.hpp>
#include <system_error>

// get_user_category

class user_category: public boost::system::error_category
{
public:

    virtual const char * name() const noexcept
    {
        return "user";
    }

    virtual std::string message( int ev ) const
    {
        char buffer[ 256 ];
        boost::core::snprintf( buffer, sizeof( buffer ), "user message %d", ev );

        return buffer;
    }
};

boost::system::error_category const & get_user_category()
{
    static user_category instance;
    return instance;
}

//

bool init_lwt = (boost::core::lwt_init(), true);

std::error_category const & cat = get_user_category();

int main()
{
    BOOST_TEST_CSTR_EQ( cat.name(), "user" );
    return boost::report_errors();
}
