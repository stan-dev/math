// Copyright 2022 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#if defined(__GNUC__) && __GNUC__ == 11
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

#include <boost/charconv.hpp>

int main()
{
    char buffer[ 32 ];
    auto r = boost::charconv::to_chars( buffer, buffer + sizeof( buffer ), 1048576 );

    int v = 0;
    boost::charconv::from_chars( buffer, r.ptr, v );
}
