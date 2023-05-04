//
// Copyright (c) 2022 Alexander Grund
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/locale/util.hpp>

#include <boost/assert.hpp>
#include <iostream>

int main()
{
    // Simple test that including a library header and using an exported function works
    std::cout << "create_utf8_converter_unique_ptr" << std::endl;
    std::unique_ptr<boost::locale::util::base_converter> cvt = boost::locale::util::create_utf8_converter_unique_ptr();
    std::cout << "create_utf8_converter" << std::endl;
    cvt = boost::locale::util::create_utf8_converter();
    cvt.reset(boost::locale::util::create_utf8_converter_new_ptr());

    if(cvt) {
        std::cout << "Created..." << std::endl;
        BOOST_ASSERT(cvt->is_thread_safe());
        BOOST_ASSERT(cvt->max_len() == 4);
    } else {
        std::cout << "Failed creation..." << std::endl;
    }

    return cvt ? 0 : 1;
}
