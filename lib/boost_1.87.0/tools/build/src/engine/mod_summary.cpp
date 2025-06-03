/*
Copyright 2024 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#include "mod_summary.h"

#include "output.h"

namespace b2 {

void summary::group(value_ref group)
{
    group_order.push_back(group);
    groups.emplace(group, group_t(new group_t::element_type));
}

void summary::message(value_ref group, value_ref message)
{
    groups[group]->push_back(message);
}

int summary::count(value_ref group) { return (int)(groups[group]->size()); }

void summary::print(value_ref group, value_ref format)
{
    std::string format_str = format;
    auto & g = groups[group];
    std::sort(g->begin(), g->end(), [](value_ref a, value_ref b) -> bool {
        return std::strcmp(a->str(), b->str()) < 0;
    });
    for (auto const & m : *g)
    {
        std::string m_str = m;
        out_printf(format->str(), m_str.c_str());
    }
}

} // namespace b2
