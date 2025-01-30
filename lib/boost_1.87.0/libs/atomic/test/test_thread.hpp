//  Copyright (c) 2023 Andrey Semashev
//
//  Distributed under the Boost Software License, Version 1.0.
//  See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_ATOMIC_TEST_THREAD_HPP_INCLUDED_
#define BOOST_ATOMIC_TEST_THREAD_HPP_INCLUDED_

#include <chrono>
#include <thread>
#include <mutex>
#include <condition_variable>
#include "test_clock.hpp"

//! Test thread class with the ability to join the thread with a timeout
class test_thread
{
private:
    std::mutex m_mutex;
    std::condition_variable m_cond;
    bool m_finished;
    std::thread m_thread;

public:
    template< typename Func >
    explicit test_thread(Func&& func) :
        m_finished(false),
        m_thread([this, func]()
        {
            try
            {
                func();
            }
            catch (...)
            {
                mark_finished();
                throw;
            }
            mark_finished();
        })
    {
    }

    test_thread(test_thread const&) = delete;
    test_thread& operator= (test_thread const&) = delete;

    void join()
    {
        m_thread.join();
    }

    template< typename Rep, typename Period >
    bool try_join_for(std::chrono::duration< Rep, Period > dur)
    {
        return try_join_until(steady_clock::now() + dur);
    }

    template< typename Clock, typename Duration >
    bool try_join_until(std::chrono::time_point< Clock, Duration > timeout)
    {
        {
            std::unique_lock< std::mutex > lock(m_mutex);
            while (!m_finished)
            {
                if (m_cond.wait_until(lock, timeout) == std::cv_status::timeout)
                    return false;
            }
        }

        join();
        return true;
    }

private:
    void mark_finished()
    {
        std::lock_guard< std::mutex > lock(m_mutex);
        m_finished = true;
        m_cond.notify_all();
    }
};

#endif // BOOST_ATOMIC_TEST_THREAD_HPP_INCLUDED_
