/*  Copyright 2019 Rene Rivera
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
 */

#ifndef B2_OPTVAL_H
#define B2_OPTVAL_H

namespace b2
{
template <typename Value>
class optval
{
public:
    typedef Value value_type;

    optval() = default;
    optval(const optval &) = default;

    template <class T>
    optval(T v) : value_is_set(true), value(v) {}

    bool has_value() const { return value_is_set; }

    operator value_type() const { return value; }

    value_type *operator->() { return &value; }

    value_type &operator*() { return value; }

    template <class T>
    optval &operator=(T v)
    {
        value = v;
        return *this;
    }

private:
    bool value_is_set = false;
    value_type value{};
};
} // namespace b2

#endif
