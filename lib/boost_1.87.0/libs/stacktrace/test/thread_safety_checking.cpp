// Copyright Antony Polukhin, 2016-2024.
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "test_impl.hpp"
#include <boost/stacktrace/stacktrace_fwd.hpp>

#include <chrono> 
#include <sstream>
#include <thread>

#include <boost/stacktrace.hpp>
#include <boost/optional.hpp>
#include <boost/core/lightweight_test.hpp>

using boost::stacktrace::stacktrace;


void main_test_loop() {
    std::size_t loops = 100;
    int Depth = 25;

    boost::optional<std::pair<stacktrace, stacktrace> > ethalon;
    std::stringstream ss_ethalon;

    while (--loops) {
        std::pair<stacktrace, stacktrace> res = function_from_library(Depth, function_from_main_translation_unit);
        if (ethalon) {
            BOOST_TEST(res == *ethalon);

            std::stringstream ss;
            ss << res.first;
            BOOST_TEST(ss.str() == ss_ethalon.str());
        } else {
            ethalon = res;
            ss_ethalon << ethalon->first;
        }
    }
}

int main() {
    const auto t = std::chrono::steady_clock::now();

    std::thread t1(main_test_loop);
    std::thread t2(main_test_loop);
    std::thread t3(main_test_loop);
    main_test_loop();

    t1.join();
    t2.join();
    t3.join();

    const auto elapsed = t - std::chrono::steady_clock::now();
    std::cout << "Elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>(
        elapsed
    ). count() << "ms";
    return boost::report_errors();
}
