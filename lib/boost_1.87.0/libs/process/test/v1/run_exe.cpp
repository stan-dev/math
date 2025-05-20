// Copyright (c) 2006, 2007 Julio M. Merino Vidal
// Copyright (c) 2008 Ilya Sokolov, Boris Schaeling
// Copyright (c) 2009 Boris Schaeling
// Copyright (c) 2010 Felipe Tanus, Boris Schaeling
// Copyright (c) 2011, 2012 Jeff Flinn, Boris Schaeling
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/core/lightweight_test.hpp>
#include <iostream>


#include <boost/process/v1/cmd.hpp>
#include <boost/process/v1/error.hpp>
#include <boost/process/v1/child.hpp>

namespace bp = boost::process;

int main(int argc, char* argv[])
{
    std::error_code ec;
    BOOST_TEST(!ec);

    bp::child c(argv[1], ec);
    BOOST_TEST(!ec);

    if (ec)
        std::cerr << ec.message() << std::endl;

    auto c2 = bp::child("doesnt-exist", ec);
    BOOST_TEST(ec);


    try
    {
        bp::child c("doesnt-exist");
        BOOST_TEST(false);
    }
    catch(bp::process_error & se)
    {
        BOOST_TEST(true);
    }

    return boost::report_errors();
}
