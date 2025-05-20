// Copyright 2018-2024 Emil Dotchevski and Reverge Studios, Inc.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/leaf/config.hpp>
#include "boost/core/lightweight_test.hpp"
#include <exception>

#ifdef BOOST_LEAF_NO_EXCEPTIONS
namespace boost
{
    [[noreturn]] void throw_exception( std::exception const & e )
    {
        std::cerr << "Terminating due to a C++ exception under BOOST_LEAF_NO_EXCEPTIONS: " << e.what();
        std::terminate();
    }

    struct source_location;
    [[noreturn]] void throw_exception( std::exception const & e, boost::source_location const & )
    {
        throw_exception(e);
    }
}
#endif
