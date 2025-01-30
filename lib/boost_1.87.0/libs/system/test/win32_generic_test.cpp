// Copyright 2022 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/system.hpp>
#include <boost/config.hpp>
#include <boost/config/pragma_message.hpp>

#if !defined(BOOST_MSSTL_VERSION) || BOOST_MSSTL_VERSION < 140

BOOST_PRAGMA_MESSAGE( "Skipping test, BOOST_MSSTL_VERSION is not defined or is less than 140" )
int main() {}

#else

#include <system_error>
#include <iostream>

int main()
{
    namespace sys = boost::system;

    int n = 0;

    for( int i = 0; i < 65000; ++i )
    {
        sys::error_code ec1( i, sys::system_category() );
        sys::error_condition en1 = ec1.default_error_condition();

        std::error_code ec2( i, std::system_category() );
        std::error_condition en2 = ec2.default_error_condition();

        if( en1 != en2 )
        {
            std::cout << i << ": " << en1 << " (" << en1.message() << ") != cond:" << en2.category().name() << ":" << en2.value() << " (" << en2.message() << ")\n";

            if( en2.category() == std::generic_category() && i != 123 ) // msvc-14.0, msvc-14.1 disagree with us on ERROR_INVALID_NAME
            {
                ++n;
            }
        }
    }

    return n < 256 ? n: 255;
}

#endif
