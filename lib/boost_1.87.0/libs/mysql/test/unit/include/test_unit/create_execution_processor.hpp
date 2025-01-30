//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CREATE_EXECUTION_PROCESSOR_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CREATE_EXECUTION_PROCESSOR_HPP

#include <boost/mysql/column_type.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/string_view.hpp>
#include <boost/mysql/throw_on_error.hpp>

#include <boost/mysql/detail/access.hpp>
#include <boost/mysql/detail/execution_processor/execution_processor.hpp>
#include <boost/mysql/detail/resultset_encoding.hpp>

#include <cstddef>

#include "test_unit/create_meta.hpp"
#include "test_unit/create_row_message.hpp"

namespace boost {
namespace mysql {
namespace test {

inline void add_meta(detail::execution_processor& proc, const std::vector<detail::coldef_view>& meta)
{
    diagnostics diag;
    proc.on_num_meta(meta.size());
    for (const auto& m : meta)
    {
        auto err = proc.on_meta(m, diag);
        throw_on_error(err, diag);
    }
}

inline void add_meta(detail::execution_processor& proc, const std::vector<column_type>& types)
{
    diagnostics diag;
    proc.on_num_meta(types.size());
    for (auto type : types)
    {
        auto err = proc.on_meta(meta_builder().type(type).build_coldef(), diag);
        throw_on_error(err, diag);
    }
}

// This is only applicable for results types (not for execution_state types)
template <class... T>
void add_row(detail::execution_processor& proc, const T&... args)
{
    auto serialized_row = create_text_row_body(args...);
    std::vector<field_view> fields;
    proc.on_row_batch_start();
    auto err = proc.on_row(serialized_row, detail::output_ref(), fields);
    throw_on_error(err);
    proc.on_row_batch_finish();
}

inline void add_ok(detail::execution_processor& proc, const detail::ok_view& pack)
{
    diagnostics diag;
    error_code err;
    if (proc.is_reading_head())
        err = proc.on_head_ok_packet(pack, diag);
    else
        err = proc.on_row_ok_packet(pack);
    throw_on_error(err, diag);
}

template <class T>
auto get_iface(T& obj) -> decltype(detail::access::get_impl(obj).get_interface())
{
    return detail::access::get_impl(obj).get_interface();
}

// Generic facility to manipulate any execution processor
// It's not owning to support any processor type without templates
class exec_access
{
    detail::execution_processor& res_;

public:
    exec_access(detail::execution_processor& obj) noexcept : res_(obj) {}

    exec_access& reset(
        detail::resultset_encoding enc = detail::resultset_encoding::text,
        metadata_mode mode = metadata_mode::minimal
    )
    {
        res_.reset(enc, mode);
        return *this;
    }
    exec_access& seqnum(std::uint8_t v)
    {
        res_.sequence_number() = v;
        return *this;
    }
    exec_access& meta(const std::vector<column_type>& types)
    {
        add_meta(res_, types);
        return *this;
    }
    exec_access& meta(const std::vector<detail::coldef_view>& meta)
    {
        add_meta(res_, meta);
        return *this;
    }
    template <class... Args>
    exec_access& row(const Args&... args)
    {
        add_row(res_, args...);
        return *this;
    }
    exec_access& ok(const detail::ok_view& pack)
    {
        add_ok(res_, pack);
        return *this;
    }
};

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
