// Copyright (c) 2022 Klemens D. Morgenstern
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/cobalt/concepts.hpp>


struct foo_aw
{
    bool await_ready();
    bool await_suspend(std::coroutine_handle<void>);
    void await_resume();
};

struct foo
{
    foo_aw operator co_await() const;
};


struct my_promise ;


struct special_aw
{
    bool await_ready();
    bool await_suspend(std::coroutine_handle<my_promise>);
    void await_resume();
};



static_assert(boost::cobalt::awaitable_type<std::suspend_never>);
static_assert(boost::cobalt::awaitable_type<foo_aw>);
static_assert(!boost::cobalt::awaitable_type<foo>);

static_assert(boost::cobalt::awaitable<std::suspend_never>);
static_assert(boost::cobalt::awaitable<foo_aw>);
static_assert(boost::cobalt::awaitable<foo>);

static_assert(boost::cobalt::awaitable<std::suspend_never, int>);
static_assert(boost::cobalt::awaitable<foo_aw, std::noop_coroutine_promise>);
static_assert(boost::cobalt::awaitable<foo, std::noop_coroutine_promise>);

static_assert(!boost::cobalt::awaitable<special_aw, std::noop_coroutine_promise>);
static_assert(!boost::cobalt::awaitable<special_aw>);
static_assert(boost::cobalt::awaitable<special_aw, my_promise>);