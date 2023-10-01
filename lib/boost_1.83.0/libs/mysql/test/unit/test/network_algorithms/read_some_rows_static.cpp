//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/detail/config.hpp>

#ifdef BOOST_MYSQL_CXX14

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/static_execution_state.hpp>

#include <boost/mysql/detail/network_algorithms.hpp>

#include <boost/mysql/impl/internal/channel/channel.hpp>

#include <boost/core/span.hpp>
#include <boost/test/unit_test.hpp>

#include "test_unit/create_channel.hpp"
#include "test_unit/create_execution_processor.hpp"
#include "test_unit/create_frame.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_ok.hpp"
#include "test_unit/create_row_message.hpp"
#include "test_unit/unit_netfun_maker.hpp"

using namespace boost::mysql::test;
using namespace boost::mysql;
using boost::span;
using boost::mysql::detail::channel;

BOOST_AUTO_TEST_SUITE(test_read_some_rows_static)

using row1 = std::tuple<int, float>;
using row2 = std::tuple<double>;

using state_t = static_execution_state<row1, row1, row2, row1, row2>;
using netfun_maker_row1 = netfun_maker_fn<std::size_t, channel&, state_t&, span<row1> >;
using netfun_maker_row2 = netfun_maker_fn<std::size_t, channel&, state_t&, span<row2> >;

struct
{
    typename netfun_maker_row1::signature read_some_rows_row1;
    typename netfun_maker_row2::signature read_some_rows_row2;
    const char* name;
} all_fns[] = {
    {netfun_maker_row1::sync_errc(&detail::read_some_rows_static_interface),
     netfun_maker_row2::sync_errc(&detail::read_some_rows_static_interface),
     "sync" },
    {netfun_maker_row1::async_errinfo(&detail::async_read_some_rows_static_interface),
     netfun_maker_row2::async_errinfo(&detail::async_read_some_rows_static_interface),
     "async"},
};

struct fixture
{
    state_t st;
    channel chan{create_channel()};
    std::array<row1, 3> storage1;
    std::array<row2, 3> storage2;

    test_stream& stream() noexcept { return get_stream(chan); }

    void add_ok() { ::add_ok(get_iface(st), ok_builder().more_results(true).build()); }

    void add_meta_row1()
    {
        add_meta(
            get_iface(st),
            {
                meta_builder().type(column_type::int_).nullable(false).build_coldef(),
                meta_builder().type(column_type::float_).nullable(false).build_coldef(),
            }
        );
    }

    void add_meta_row2()
    {
        add_meta(
            get_iface(st),
            {
                meta_builder().type(column_type::double_).nullable(false).build_coldef(),
            }
        );
    }
};

BOOST_AUTO_TEST_CASE(repeated_row_types)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.add_meta_row1();

            // 1st resultset: row1
            fix.stream()
                .add_bytes(create_text_row_message(0, 10, 4.2f))
                .add_bytes(create_text_row_message(1, 11, 4.3f));

            std::size_t num_rows = fns.read_some_rows_row1(fix.chan, fix.st, fix.storage1).get();
            BOOST_TEST_REQUIRE(num_rows == 2u);
            BOOST_TEST((fix.storage1[0] == row1{10, 4.2f}));
            BOOST_TEST((fix.storage1[1] == row1{11, 4.3f}));

            // Advance resultset
            fix.add_ok();
            fix.add_meta_row1();
            BOOST_TEST_REQUIRE(fix.st.should_read_rows());

            // 2nd resultset: row1 again
            fix.stream().add_bytes(create_text_row_message(2, 13, 0.2f));
            num_rows = fns.read_some_rows_row1(fix.chan, fix.st, fix.storage1).get();
            BOOST_TEST_REQUIRE(num_rows == 1u);
            BOOST_TEST((fix.storage1[0] == row1{13, 0.2f}));

            // Advance resultset
            fix.add_ok();
            fix.add_meta_row2();
            BOOST_TEST_REQUIRE(fix.st.should_read_rows());

            // 3rd resultset: row2
            fix.stream().add_bytes(create_text_row_message(3, 9.1));
            num_rows = fns.read_some_rows_row2(fix.chan, fix.st, fix.storage2).get();
            BOOST_TEST_REQUIRE(num_rows == 1u);
            BOOST_TEST((fix.storage2[0] == row2{9.1}));

            // Advance resultset
            fix.add_ok();
            fix.add_meta_row1();
            BOOST_TEST_REQUIRE(fix.st.should_read_rows());

            // 4th resultset: row1
            fix.stream().add_bytes(create_text_row_message(4, 43, 0.7f));
            num_rows = fns.read_some_rows_row1(fix.chan, fix.st, fix.storage1).get();
            BOOST_TEST_REQUIRE(num_rows == 1u);
            BOOST_TEST((fix.storage1[0] == row1{43, 0.7f}));

            // Advance resultset
            fix.add_ok();
            fix.add_meta_row2();
            BOOST_TEST_REQUIRE(fix.st.should_read_rows());

            // 5th resultset: row2
            fix.stream().add_bytes(create_text_row_message(5, 99.9));
            num_rows = fns.read_some_rows_row2(fix.chan, fix.st, fix.storage2).get();
            BOOST_TEST_REQUIRE(num_rows == 1u);
            BOOST_TEST((fix.storage2[0] == row2{99.9}));
        }
    }
}

BOOST_AUTO_TEST_CASE(error_row_type_mismatch)
{
    for (const auto& fns : all_fns)
    {
        BOOST_TEST_CONTEXT(fns.name)
        {
            fixture fix;
            fix.add_meta_row1();

            // 1st resultset: row1. Note that this will consume the message
            fix.stream().add_bytes(create_text_row_message(0, 10, 4.2f));
            fns.read_some_rows_row2(fix.chan, fix.st, fix.storage2)
                .validate_error_exact(client_errc::row_type_mismatch);

            // Advance resultset
            fix.add_ok();
            fix.add_meta_row1();
            fix.add_ok();
            fix.add_meta_row2();
            BOOST_TEST_REQUIRE(fix.st.should_read_rows());

            // 3rd resultset: row2
            fix.stream().add_bytes(create_text_row_message(1, 9.1));
            fns.read_some_rows_row1(fix.chan, fix.st, fix.storage1)
                .validate_error_exact(client_errc::row_type_mismatch);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

#endif
