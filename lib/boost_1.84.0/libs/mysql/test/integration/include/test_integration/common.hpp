//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_COMMON_HPP
#define BOOST_MYSQL_TEST_INTEGRATION_INCLUDE_TEST_INTEGRATION_COMMON_HPP

#include <boost/mysql/execution_state.hpp>
#include <boost/mysql/handshake_params.hpp>
#include <boost/mysql/metadata_collection_view.hpp>
#include <boost/mysql/results.hpp>

#include <boost/asio/io_context.hpp>
#include <boost/asio/ssl/context.hpp>
#include <boost/asio/system_executor.hpp>
#include <boost/test/unit_test.hpp>

#include <initializer_list>
#include <thread>

#include "test_integration/er_connection.hpp"
#include "test_integration/er_network_variant.hpp"
#include "test_integration/metadata_validator.hpp"
#include "test_integration/network_test.hpp"

namespace boost {
namespace mysql {
namespace test {

struct network_fixture_base
{
    handshake_params params{"integ_user", "integ_password", "boost_mysql_integtests"};
    boost::asio::io_context ctx;
    boost::asio::ssl::context ssl_ctx{boost::asio::ssl::context::tls_client};
};

struct network_fixture : network_fixture_base
{
    er_network_variant* var{};
    er_connection_ptr conn;

    ~network_fixture()
    {
        if (conn)
        {
            conn->sync_close();
        }
    }

    void setup(er_network_variant* variant)
    {
        var = variant;
        conn = var->create_connection(ctx.get_executor(), ssl_ctx);
        conn->set_metadata_mode(metadata_mode::full);
    }

    void setup_and_physical_connect(er_network_variant* net)
    {
        setup(net);
        conn->physical_connect();
    }

    void setup_and_connect(er_network_variant* variant, ssl_mode m = ssl_mode::require)
    {
        setup(variant);
        connect(m);
    }

    void set_credentials(string_view user, string_view password)
    {
        params.set_username(user);
        params.set_password(password);
    }

    // Verifies that we are or are not using SSL, depending on whether the stream supports it or not
    void validate_ssl(ssl_mode m = ssl_mode::require)
    {
        bool expected = (m == ssl_mode::require || m == ssl_mode::enable) && var->supports_ssl();
        BOOST_TEST(conn->uses_ssl() == expected);
    }

    void handshake(ssl_mode m = ssl_mode::require)
    {
        assert(conn);
        params.set_ssl(m);
        conn->handshake(params).validate_no_error();
        validate_ssl(m);
    }

    void connect(ssl_mode m = ssl_mode::require)
    {
        assert(conn);
        params.set_ssl(m);
        conn->connect(params).validate_no_error();
        validate_ssl(m);
    }

    // Call this in the fixture setup of any test invoking write
    // operations on the database, to prevent race conditions,
    // make the testing environment more stable and speed up the tests
    void start_transaction()
    {
        results result;
        conn->execute("START TRANSACTION", result).get();
    }
};

// To be used as sample in data driven tests, when a test case should be run
// over all different network_function's.
struct network_sample
{
    er_network_variant* net;

    network_sample(er_network_variant* var) : net(var) {}

    void set_test_attributes(boost::unit_test::test_case& test) const
    {
        if (net->supports_ssl())
        {
            test.add_label("ssl");
        }
        if (net->is_unix_socket())
        {
            test.add_label("unix");
        }
        test.add_label(net->variant_name());
    }
};

inline std::ostream& operator<<(std::ostream& os, const network_sample& value)
{
    return os << value.net->stream_name() << "_" << value.net->variant_name();
}

inline std::vector<network_sample> create_network_samples(std::initializer_list<const char*> names)
{
    std::vector<network_sample> res;
    for (const char* name : names)
        res.emplace_back(get_variant(name));
    return res;
}

inline std::vector<network_sample> create_all_network_samples()
{
    std::vector<network_sample> res;
    for (auto* var : all_variants())
        res.emplace_back(var);
    return res;
}

inline const std::vector<network_sample>& all_network_samples()
{
    static std::vector<network_sample> res = create_all_network_samples();
    return res;
}

inline void validate_2fields_meta(const metadata_collection_view& fields, const std::string& table)

{
    validate_meta(
        fields,
        {meta_validator(table, "id", column_type::int_),
         meta_validator(table, "field_varchar", column_type::varchar)}
    );
}

inline void validate_eof(
    const execution_state& st,
    unsigned affected_rows = 0,
    unsigned warnings = 0,
    unsigned last_insert = 0,
    string_view info = ""
)
{
    BOOST_TEST_REQUIRE(st.complete());
    BOOST_TEST(st.affected_rows() == affected_rows);
    BOOST_TEST(st.warning_count() == warnings);
    BOOST_TEST(st.last_insert_id() == last_insert);
    BOOST_TEST(st.info() == info);
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif /* TEST_INTEGRATION_INTEGRATION_TEST_COMMON_HPP_ */
