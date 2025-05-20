// Copyright 2018-2024 Emil Dotchevski and Reverge Studios, Inc.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// This is a simple program that demonstrates the use of LEAF to transport error
// objects between threads, without using exception handling. See try_capture_all_exceptions.cpp
// for the version that uses exception handling.

#include <boost/leaf/config.hpp>

#if !BOOST_LEAF_CFG_CAPTURE

#include <iostream>

int main()
{
    std::cout << "Unit test not applicable." << std::endl;
    return 0;
}

#else

#include <boost/leaf.hpp>
#include <vector>
#include <string>
#include <future>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <thread>

namespace leaf = boost::leaf;

// Define several error types.
struct e_thread_id { std::thread::id value; };
struct e_failure_info1 { std::string value; };
struct e_failure_info2 { int value; };

// A type that represents a successfully returned result from a task.
struct task_result { };

 // This is our task function. It produces objects of type task_result, but it
 // may fail.
leaf::result<task_result> task()
{
    bool succeed = (rand()%4) != 0; //...at random.
    if( succeed )
        return { };
    else
        return leaf::new_error(
            e_thread_id{std::this_thread::get_id()},
            e_failure_info1{"info"},
            e_failure_info2{42} );
};

int main()
{
    int const task_count = 42;

    // Container to collect the generated std::future objects.
    std::vector<std::future<leaf::result<task_result>>> fut;

    // Launch the tasks, but rather than launching the task function directly,
    // we use leaf::try_capture_all: in case of a failure, the returned leaf::result<>
    // will capture all error objects.
    std::generate_n( std::back_inserter(fut), task_count,
        [&]
        {
            return std::async(
                std::launch::async,
                [&]
                {
                    return leaf::try_capture_all(task);
                } );
        } );

    // Wait on the futures, get the task results, handle errors.
    for( auto & f : fut )
    {
        f.wait();

        leaf::try_handle_all(
            [&]() -> leaf::result<void>
            {
                BOOST_LEAF_AUTO(r, f.get());

                // Success! Use r to access task_result.
                std::cout << "Success!" << std::endl;
                (void) r; // Presumably we'll somehow use the task_result.
                return { };
            },
            []( e_failure_info1 const & v1, e_failure_info2 const & v2, e_thread_id const & tid )
            {
                std::cerr << "Error in thread " << tid.value << "! failure_info1: " << v1.value << ", failure_info2: " << v2.value << std::endl;
            },
            []( leaf::diagnostic_info const & unmatched )
            {
                std::cerr <<
                    "Unknown failure detected" << std::endl <<
                    "Cryptic diagnostic information follows" << std::endl <<
                    unmatched;
            } );
    }
}

////////////////////////////////////////

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

#endif
