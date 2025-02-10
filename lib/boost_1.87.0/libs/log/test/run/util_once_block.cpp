/*
 *          Copyright Andrey Semashev 2007 - 2015.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 */
/*!
 * \file   util_once_block.cpp
 * \author Andrey Semashev
 * \date   24.06.2010
 *
 * \brief  This header contains tests for once-blocks implementation.
 *
 * The test was adopted from test_once.cpp test of Boost.Thread.
 */

#define BOOST_TEST_MODULE util_once_block

#include <boost/log/utility/once_block.hpp>
#include <boost/test/unit_test.hpp>

#if !defined(BOOST_LOG_NO_THREADS)

#include <mutex>
#include <thread>
#include "test_barrier.hpp"

namespace logging = boost::log;

enum config
{
    THREAD_COUNT = 20,
    LOOP_COUNT = 100
};

std::mutex m;
typedef std::lock_guard< std::mutex > scoped_lock;

logging::once_block_flag flag = BOOST_LOG_ONCE_BLOCK_INIT;
int var_to_init_once_flag = 0;

void initialize_variable()
{
    // ensure that if multiple threads get in here, they are serialized, so we can see the effect
    scoped_lock lock(m);
    ++var_to_init_once_flag;
}


void once_block_flag_thread(test_barrier& barrier)
{
    int my_once_value = 0;
    barrier.arrive_and_wait();
    for (unsigned int i = 0; i < LOOP_COUNT; ++i)
    {
        BOOST_LOG_ONCE_BLOCK_FLAG(flag)
        {
            initialize_variable();
        }

        my_once_value = var_to_init_once_flag;
        if (my_once_value != 1)
        {
            break;
        }
    }
    scoped_lock lock(m);
    BOOST_CHECK_EQUAL(my_once_value, 1);
}

// The test checks if the BOOST_LOG_ONCE_BLOCK_FLAG macro works
BOOST_AUTO_TEST_CASE(once_block_flag)
{
    std::thread threads[THREAD_COUNT];
    test_barrier barrier(static_cast< unsigned int >(THREAD_COUNT));

    try
    {
        for (unsigned int i = 0; i < THREAD_COUNT; ++i)
            threads[i] = std::thread([&barrier]() { once_block_flag_thread(barrier); });

        for (unsigned int i = 0; i < THREAD_COUNT; ++i)
            threads[i].join();
    }
    catch (...)
    {
        barrier.wake_all();
        for (unsigned int i = 0; i < THREAD_COUNT; ++i)
        {
            if (threads[i].joinable())
                threads[i].join();
        }
        throw;
    }

    BOOST_CHECK_EQUAL(var_to_init_once_flag, 1);
}

int var_to_init_once = 0;

void once_block_thread(test_barrier& barrier)
{
    int my_once_value = 0;
    barrier.arrive_and_wait();
    for (unsigned int i = 0; i < LOOP_COUNT; ++i)
    {
        BOOST_LOG_ONCE_BLOCK()
        {
            scoped_lock lock(m);
            ++var_to_init_once;
        }

        my_once_value = var_to_init_once;
        if (my_once_value != 1)
        {
            break;
        }
    }

    scoped_lock lock(m);
    BOOST_CHECK_EQUAL(my_once_value, 1);
}

// The test checks if the BOOST_LOG_ONCE_BLOCK macro works
BOOST_AUTO_TEST_CASE(once_block)
{
    std::thread threads[THREAD_COUNT];
    test_barrier barrier(static_cast< unsigned int >(THREAD_COUNT));

    try
    {
        for (unsigned int i = 0; i < THREAD_COUNT; ++i)
            threads[i] = std::thread([&barrier]() { once_block_thread(barrier); });

        for (unsigned int i = 0; i < THREAD_COUNT; ++i)
            threads[i].join();
    }
    catch (...)
    {
        barrier.wake_all();
        for (unsigned int i = 0; i < THREAD_COUNT; ++i)
        {
            if (threads[i].joinable())
                threads[i].join();
        }
        throw;
    }

    BOOST_CHECK_EQUAL(var_to_init_once, 1);
}


struct my_exception
{
};

unsigned int pass_counter = 0;
unsigned int exception_counter = 0;

void once_block_with_exception_thread(test_barrier& barrier)
{
    barrier.arrive_and_wait();
    try
    {
        BOOST_LOG_ONCE_BLOCK()
        {
            scoped_lock lock(m);
            ++pass_counter;
            if (pass_counter < 3)
            {
                throw my_exception();
            }
        }
    }
    catch (my_exception&)
    {
        scoped_lock lock(m);
        ++exception_counter;
    }
}

// The test verifies that the once_block flag is not set if an exception is thrown from the once-block
BOOST_AUTO_TEST_CASE(once_block_retried_on_exception)
{
    std::thread threads[THREAD_COUNT];
    test_barrier barrier(static_cast< unsigned int >(THREAD_COUNT));

    try
    {
        for (unsigned int i = 0; i < THREAD_COUNT; ++i)
            threads[i] = std::thread([&barrier]() { once_block_with_exception_thread(barrier); });

        for (unsigned int i = 0; i < THREAD_COUNT; ++i)
            threads[i].join();
    }
    catch (...)
    {
        barrier.wake_all();
        for (unsigned int i = 0; i < THREAD_COUNT; ++i)
        {
            if (threads[i].joinable())
                threads[i].join();
        }
        throw;
    }

    BOOST_CHECK_EQUAL(pass_counter, 3u);
    BOOST_CHECK_EQUAL(exception_counter, 2u);
}

#else // BOOST_LOG_NO_THREADS

namespace logging = boost::log;

enum config
{
    LOOP_COUNT = 100
};

logging::once_block_flag flag = BOOST_LOG_ONCE_BLOCK_INIT;
int var_to_init_once_flag = 0;

void initialize_variable()
{
    ++var_to_init_once_flag;
}


// The test checks if the BOOST_LOG_ONCE_BLOCK_FLAG macro works
BOOST_AUTO_TEST_CASE(once_block_flag)
{
    for (unsigned int i = 0; i < LOOP_COUNT; ++i)
    {
        BOOST_LOG_ONCE_BLOCK_FLAG(flag)
        {
            initialize_variable();
        }

        if (var_to_init_once_flag != 1)
        {
            break;
        }
    }

    BOOST_CHECK_EQUAL(var_to_init_once_flag, 1);
}

int var_to_init_once = 0;

// The test checks if the BOOST_LOG_ONCE_BLOCK macro works
BOOST_AUTO_TEST_CASE(once_block)
{
    for (unsigned int i = 0; i < LOOP_COUNT; ++i)
    {
        BOOST_LOG_ONCE_BLOCK()
        {
            ++var_to_init_once;
        }

        if (var_to_init_once != 1)
        {
            break;
        }
    }

    BOOST_CHECK_EQUAL(var_to_init_once, 1);
}

struct my_exception
{
};

unsigned int pass_counter = 0;
unsigned int exception_counter = 0;

// The test verifies that the once_block flag is not set if an exception is thrown from the once-block
BOOST_AUTO_TEST_CASE(once_block_retried_on_exception)
{
    for (unsigned int i = 0; i < LOOP_COUNT; ++i)
    {
        try
        {
            BOOST_LOG_ONCE_BLOCK()
            {
                ++pass_counter;
                if (pass_counter < 3)
                {
                    throw my_exception();
                }
            }
        }
        catch (my_exception&)
        {
            ++exception_counter;
        }
    }

    BOOST_CHECK_EQUAL(pass_counter, 3u);
    BOOST_CHECK_EQUAL(exception_counter, 2u);
}

#endif // BOOST_LOG_NO_THREADS
