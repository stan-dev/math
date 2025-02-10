//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/column_type.hpp>
#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/metadata.hpp>
#include <boost/mysql/metadata_mode.hpp>
#include <boost/mysql/string_view.hpp>
#include <boost/mysql/throw_on_error.hpp>

#include <boost/mysql/detail/coldef_view.hpp>
#include <boost/mysql/detail/execution_processor/execution_processor.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/check_meta.hpp"
#include "test_common/printing.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/mock_execution_processor.hpp"
#include "test_unit/printing.hpp"

using namespace boost::mysql::detail;
using namespace boost::mysql::test;
using boost::mysql::column_type;
using boost::mysql::diagnostics;
using boost::mysql::error_code;
using boost::mysql::metadata_mode;
using boost::mysql::string_view;
using boost::mysql::throw_on_error;

namespace {

class spy_execution_processor : public mock_execution_processor
{
public:
    struct
    {
        bool is_last;
    } on_meta_call;

private:
    error_code on_meta_impl(const coldef_view& coldef, bool is_last, diagnostics& diag) override
    {
        on_meta_call.is_last = is_last;
        return mock_execution_processor::on_meta_impl(coldef, is_last, diag);
    }
};

void check_reading_first(const execution_processor& st)
{
    BOOST_TEST(st.is_reading_first());
    BOOST_TEST(!st.is_reading_first_subseq());
    BOOST_TEST(st.is_reading_head());
    BOOST_TEST(!st.is_reading_meta());
    BOOST_TEST(!st.is_reading_rows());
    BOOST_TEST(!st.is_complete());
}

void check_reading_first_subseq(const execution_processor& st)
{
    BOOST_TEST(!st.is_reading_first());
    BOOST_TEST(st.is_reading_first_subseq());
    BOOST_TEST(st.is_reading_head());
    BOOST_TEST(!st.is_reading_meta());
    BOOST_TEST(!st.is_reading_rows());
    BOOST_TEST(!st.is_complete());
}

void check_reading_meta(const execution_processor& st)
{
    BOOST_TEST(!st.is_reading_first());
    BOOST_TEST(!st.is_reading_first_subseq());
    BOOST_TEST(!st.is_reading_head());
    BOOST_TEST(st.is_reading_meta());
    BOOST_TEST(!st.is_reading_rows());
    BOOST_TEST(!st.is_complete());
}

void check_reading_rows(const execution_processor& st)
{
    BOOST_TEST(!st.is_reading_first());
    BOOST_TEST(!st.is_reading_first_subseq());
    BOOST_TEST(!st.is_reading_head());
    BOOST_TEST(!st.is_reading_meta());
    BOOST_TEST(st.is_reading_rows());
    BOOST_TEST(!st.is_complete());
}

void check_complete(const execution_processor& st)
{
    BOOST_TEST(!st.is_reading_first());
    BOOST_TEST(!st.is_reading_first_subseq());
    BOOST_TEST(!st.is_reading_head());
    BOOST_TEST(!st.is_reading_meta());
    BOOST_TEST(!st.is_reading_rows());
    BOOST_TEST(st.is_complete());
}

BOOST_AUTO_TEST_SUITE(test_execution_processor)

BOOST_AUTO_TEST_CASE(default_ctor)
{
    mock_execution_processor p;
    check_reading_first(p);
    BOOST_TEST(p.encoding() == resultset_encoding::text);
    BOOST_TEST(p.sequence_number() == 0u);
    BOOST_TEST(p.meta_mode() == metadata_mode::minimal);
}

BOOST_AUTO_TEST_CASE(reset)
{
    mock_execution_processor p;
    p.on_num_meta(42);
    p.sequence_number() = 42u;
    p.reset(resultset_encoding::binary, metadata_mode::full);
    check_reading_first(p);
    BOOST_TEST(p.encoding() == resultset_encoding::binary);
    BOOST_TEST(p.sequence_number() == 0u);
    BOOST_TEST(p.meta_mode() == metadata_mode::full);
}

BOOST_AUTO_TEST_CASE(states)
{
    mock_execution_processor p;
    diagnostics diag;

    check_reading_first(p);

    p.on_num_meta(1);
    check_reading_meta(p);

    auto err = p.on_meta(meta_builder().build_coldef(), diag);
    throw_on_error(err, diag);
    check_reading_rows(p);

    err = p.on_row_ok_packet(ok_builder().more_results(true).build());
    throw_on_error(err, diag);
    check_reading_first_subseq(p);

    err = p.on_head_ok_packet(ok_builder().build(), diag);
    throw_on_error(err, diag);
    check_complete(p);
}

BOOST_AUTO_TEST_CASE(on_meta_one_column)
{
    spy_execution_processor p;
    diagnostics diag;
    p.reset(resultset_encoding::text, metadata_mode::minimal);
    p.on_num_meta(1);

    auto err = p.on_meta(meta_builder().type(column_type::bit).build_coldef(), diag);

    // Verify
    BOOST_TEST(err == error_code());
    p.num_calls().reset(1).on_num_meta(1).on_meta(1).validate();
    check_meta(p.meta(), {column_type::bit});
    BOOST_TEST(p.on_meta_call.is_last);
}

BOOST_AUTO_TEST_CASE(on_meta_several_columns)
{
    spy_execution_processor p;
    diagnostics diag;
    p.reset(resultset_encoding::text, metadata_mode::full);
    p.on_num_meta(2);

    // 1st field
    auto err = p.on_meta(meta_builder().type(column_type::bit).build_coldef(), diag);

    // Verify
    BOOST_TEST(err == error_code());
    p.num_calls().reset(1).on_num_meta(1).on_meta(1).validate();
    check_meta(p.meta(), {column_type::bit});
    BOOST_TEST(!p.on_meta_call.is_last);

    // 2nd column
    err = p.on_meta(meta_builder().type(column_type::varchar).build_coldef(), diag);

    // Verify
    BOOST_TEST(err == error_code());
    p.num_calls().reset(1).on_num_meta(1).on_meta(2).validate();
    check_meta(p.meta(), {column_type::bit, column_type::varchar});
    BOOST_TEST(p.on_meta_call.is_last);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace