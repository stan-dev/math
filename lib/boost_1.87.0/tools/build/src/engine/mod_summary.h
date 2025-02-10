/*
Copyright 2024 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#ifndef B2_MOD_SUMMARY_H
#define B2_MOD_SUMMARY_H

#include "config.h"

#include "bind.h"
#include "value.h"

#include <memory>
#include <unordered_map>
#include <vector>

namespace b2 {

class summary : public object
{
    public:
    void group(value_ref group);
    void message(value_ref group, value_ref message);
    int count(value_ref group);
    void print(value_ref group, value_ref format);

    private:
    using group_t = std::unique_ptr<std::vector<value_ref>>;
    using groups_t = std::unordered_map<value_ref,
        group_t,
        value_ref::hash_function,
        value_ref::equal_function>;

    groups_t groups;
    std::vector<value_ref> group_order;
};

struct summary_module : b2::bind::module_<summary_module>
{
    const char * module_name = "summary";

    template <class Binder>
    void def(Binder & binder)
    {
        binder.def_class("summary", type_<summary>())
            .def(init_<>())
            .def(&summary::group, "group", "group" * _1)
            .def(&summary::message, "message", "group" * _1,
                "message" * _1n)
            .def(&summary::count, "count", "group" * _1)
            .def(&summary::print, "print", "group" * _1,
                "format" * _1);
    }
};

} // namespace b2

#endif
