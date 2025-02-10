// Copyright 2024 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt)

#if !defined(_WIN32)

#include <boost/config/pragma_message.hpp>

BOOST_PRAGMA_MESSAGE( "Test skipped because _WIN32 is not defined" )

#else

#include <windows.h>
#include <boost/crc.hpp>

#endif

int main()
{
}
