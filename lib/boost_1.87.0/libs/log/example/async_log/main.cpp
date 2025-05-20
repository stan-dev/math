/*
 *          Copyright Andrey Semashev 2007 - 2015.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 */
/*!
 * \file   main.cpp
 * \author Andrey Semashev
 * \date   30.08.2009
 *
 * \brief  An example of asynchronous logging in multiple threads.
 */

// #define BOOST_LOG_DYN_LINK 1

#include <stdexcept>
#include <thread>
#include <string>
#include <iostream>
#include <fstream>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/compat/latch.hpp>

#include <boost/log/common.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/attributes.hpp>
#include <boost/log/sinks.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/log/utility/record_ordering.hpp>

namespace logging = boost::log;
namespace attrs = boost::log::attributes;
namespace src = boost::log::sources;
namespace sinks = boost::log::sinks;
namespace expr = boost::log::expressions;
namespace keywords = boost::log::keywords;

using boost::shared_ptr;

enum
{
    LOG_RECORDS_TO_WRITE = 10000,
    THREAD_COUNT = 2
};

BOOST_LOG_INLINE_GLOBAL_LOGGER_DEFAULT(test_lg, src::logger_mt)

//! This function is executed in multiple threads
void thread_fun(boost::compat::latch& latch)
{
    // Wait until all threads are created
    latch.arrive_and_wait();

    // Here we go. First, identify the thread.
    BOOST_LOG_SCOPED_THREAD_TAG("ThreadID", std::this_thread::get_id());

    // Now, do some logging
    for (unsigned int i = 0; i < LOG_RECORDS_TO_WRITE; ++i)
    {
        BOOST_LOG(test_lg::get()) << "Log record " << i;
    }
}

int main(int argc, char* argv[])
{
    try
    {
        // Open a rotating text file
        shared_ptr< std::ostream > strm(new std::ofstream("test.log"));
        if (!strm->good())
            throw std::runtime_error("Failed to open a text log file");

        // Create a text file sink
        typedef sinks::text_ostream_backend backend_t;
        typedef sinks::asynchronous_sink<
            backend_t,
            sinks::unbounded_ordering_queue<
                logging::attribute_value_ordering< unsigned int, std::less< unsigned int > >
            >
        > sink_t;
        shared_ptr< sink_t > sink(new sink_t(
            boost::make_shared< backend_t >(),
            // We'll apply record ordering to ensure that records from different threads go sequentially in the file
            keywords::order = logging::make_attr_ordering< unsigned int >("RecordID", std::less< unsigned int >())));

        sink->locked_backend()->add_stream(strm);

        sink->set_formatter
        (
            expr::format("%1%: [%2%] [%3%] - %4%")
                % expr::attr< unsigned int >("RecordID")
                % expr::attr< boost::posix_time::ptime >("TimeStamp")
                % expr::attr< std::thread::id >("ThreadID")
                % expr::smessage
        );

        // Add it to the core
        logging::core::get()->add_sink(sink);

        // Add some attributes too
        logging::core::get()->add_global_attribute("TimeStamp", attrs::local_clock());
        logging::core::get()->add_global_attribute("RecordID", attrs::counter< unsigned int >());

        // Create logging threads
        boost::compat::latch latch(THREAD_COUNT);
        std::thread threads[THREAD_COUNT];
        for (unsigned int i = 0; i < THREAD_COUNT; ++i)
            threads[i] = std::thread([&latch]() { thread_fun(latch); });

        // Wait until all action ends
        for (unsigned int i = 0; i < THREAD_COUNT; ++i)
            threads[i].join();

        // Flush all buffered records
        sink->stop();
        sink->flush();

        return 0;
    }
    catch (std::exception& e)
    {
        std::cout << "FAILURE: " << e.what() << std::endl;
        return 1;
    }
}
