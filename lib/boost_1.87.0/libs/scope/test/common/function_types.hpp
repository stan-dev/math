/*
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * https://www.boost.org/LICENSE_1_0.txt)
 *
 * Copyright (c) 2023 Andrey Semashev
 */
/*!
 * \file   function_types.hpp
 * \author Andrey Semashev
 *
 * \brief  This file contains various function types for scope guard tests.
 */

#ifndef BOOST_SCOPE_TESTS_COMMON_FUNCTION_TYPES_H_INCLUDED_
#define BOOST_SCOPE_TESTS_COMMON_FUNCTION_TYPES_H_INCLUDED_

#include <stdexcept>

struct normal_func
{
    int* m_n;

    explicit normal_func(int& n) noexcept :
        m_n(&n)
    {
    }

    void operator()() const noexcept
    {
        ++(*m_n);
    }
};

struct moveable_only_func
{
    int* m_n;

    explicit moveable_only_func(int& n) noexcept :
        m_n(&n)
    {
    }

    moveable_only_func(moveable_only_func const&) = delete;
    moveable_only_func& operator= (moveable_only_func const&) = delete;

    moveable_only_func(moveable_only_func&&) = default;
    moveable_only_func& operator= (moveable_only_func&&) = default;

    void operator()() const noexcept
    {
        ++(*m_n);
    }
};

struct throw_on_copy_func
{
    int* m_n;

    explicit throw_on_copy_func(int& n) noexcept :
        m_n(&n)
    {
    }

    throw_on_copy_func(throw_on_copy_func const&)
    {
        throw std::runtime_error("throw_on_copy_func copy ctor");
    }

    throw_on_copy_func& operator= (throw_on_copy_func const&) = delete;

    void operator()() const noexcept
    {
        ++(*m_n);
    }
};

struct throw_on_move_func
{
    int* m_n;

    explicit throw_on_move_func(int& n) noexcept :
        m_n(&n)
    {
    }

    throw_on_move_func(throw_on_move_func const&) = default;
    throw_on_move_func& operator= (throw_on_move_func const&) = delete;

    throw_on_move_func(throw_on_move_func&&)
    {
        throw std::runtime_error("throw_on_move_func move ctor");
    }

    throw_on_move_func& operator= (throw_on_move_func&&) = delete;

    void operator()() const noexcept
    {
        ++(*m_n);
    }
};

struct throw_on_call_func
{
    int* m_n;
    bool* m_destroyed;

    explicit throw_on_call_func(int& n, bool& destroyed) noexcept :
        m_n(&n),
        m_destroyed(&destroyed)
    {
    }

    ~throw_on_call_func()
    {
        *m_destroyed = true;
    }

    void operator()() const
    {
        ++(*m_n);
        throw std::runtime_error("throw_on_call_func call");
    }
};

#endif // BOOST_SCOPE_TESTS_COMMON_FUNCTION_TYPES_H_INCLUDED_
