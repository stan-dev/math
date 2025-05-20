//  Copyright (c) 2024 Andrey Semashev
//
//  Distributed under the Boost Software License, Version 1.0.
//  See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_LOG_TEST_BARRIER_HPP_INCLUDED_
#define BOOST_LOG_TEST_BARRIER_HPP_INCLUDED_

#include <mutex>
#include <condition_variable>

//! A simplified version of thread barrier from Boost.Thread and C++20 std::barrier
class test_barrier
{
private:
    std::mutex m_mutex;
    std::condition_variable m_cond;
    unsigned int m_generation;
    unsigned int m_count;
    const unsigned int m_initial_count;

public:
    explicit test_barrier(unsigned int initial_count) :
        m_generation(0u), m_count(initial_count), m_initial_count(initial_count)
    {
    }

    test_barrier(test_barrier const&) = delete;
    test_barrier& operator= (test_barrier const&) = delete;

    void arrive_and_wait()
    {
        std::unique_lock< std::mutex > lock(m_mutex);

        --m_count;
        if (m_count == 0u)
        {
            ++m_generation;
            m_count = m_initial_count;
            m_cond.notify_all();
            return;
        }

        const unsigned int generation = m_generation;
        do
        {
            m_cond.wait(lock);
        }
        while (m_generation == generation);
    }

    void wake_all()
    {
        std::lock_guard< std::mutex > lock(m_mutex);
        ++m_generation;
        m_count = m_initial_count;
        m_cond.notify_all();
    }
};

#endif // BOOST_LOG_TEST_BARRIER_HPP_INCLUDED_
