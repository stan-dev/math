//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/connect_params.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/ssl_mode.hpp>
#include <boost/mysql/string_view.hpp>
#include <boost/mysql/tcp_ssl.hpp>
#include <boost/mysql/unix.hpp>
#include <boost/mysql/unix_ssl.hpp>

#include <boost/asio/ip/tcp.hpp>
#include <boost/asio/local/stream_protocol.hpp>
#include <boost/asio/ssl/context.hpp>
#include <boost/asio/ssl/host_name_verification.hpp>
#include <boost/asio/ssl/verify_mode.hpp>
#include <boost/assert/source_location.hpp>
#include <boost/test/data/monomorphic/collection.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "test_common/ci_server.hpp"
#include "test_common/create_basic.hpp"
#include "test_common/io_context_fixture.hpp"
#include "test_common/netfun_maker.hpp"
#include "test_common/network_result.hpp"
#include "test_common/source_location.hpp"
#include "test_integration/any_connection_fixture.hpp"
#include "test_integration/connect_params_builder.hpp"
#include "test_integration/server_ca.hpp"
#include "test_integration/server_features.hpp"
#include "test_integration/tcp_connection_fixture.hpp"

using namespace boost::mysql;
using namespace boost::mysql::test;
using boost::test_tools::per_element;
namespace asio = boost::asio;
namespace data = boost::unit_test::data;

namespace {

BOOST_AUTO_TEST_SUITE(test_handshake)

// Handshake is the most convoluted part of MySQL protocol,
// and is in active development in current MySQL versions.
// We try to test all combinations of auth methods/transports.
// Note that fixtures take care of closing the connection can be closed successfully
struct transport_test_case
{
    string_view name;
    connect_params params;
    bool expect_ssl;
};
std::ostream& operator<<(std::ostream& os, const transport_test_case& tc) { return os << tc.name; }

std::vector<transport_test_case> make_secure_transports()
{
    std::vector<transport_test_case> res{
        {"tcp_ssl", connect_params_builder().ssl(ssl_mode::require).build(), true},
    };

#ifdef BOOST_ASIO_HAS_LOCAL_SOCKETS
    if (get_server_features().unix_sockets)
    {
        res.push_back({"unix", connect_params_builder().set_unix().build(), false});
    }
#endif

    return res;
}

std::vector<transport_test_case> make_all_transports()
{
    auto res = make_secure_transports();
    res.push_back({"tcp", connect_params_builder().ssl(ssl_mode::disable).build(), false});
    return res;
}

auto secure_transports = make_secure_transports();
auto all_transports = make_all_transports();

// Check whether the connection is using SSL or not
template <class Conn>
void check_ssl(Conn& conn, bool expected, boost::source_location loc = BOOST_MYSQL_CURRENT_LOCATION)
{
    BOOST_TEST_CONTEXT("Called from " << loc)
    {
        // Check that the client thinks it's using SSL
        BOOST_TEST(conn.uses_ssl() == expected);

        // Check that the server is using SSL
        results r;
        conn.async_execute("SHOW STATUS LIKE 'ssl_version'", r, as_netresult).validate_no_error(loc);
        bool server_tls = r.rows().at(0).at(1).as_string().starts_with("TLS");
        BOOST_TEST(server_tls == expected);
    }
}

// mysql_native_password
BOOST_AUTO_TEST_SUITE(mysql_native_password)

constexpr const char* regular_user = "mysqlnp_user";
constexpr const char* regular_passwd = "mysqlnp_password";
constexpr const char* empty_user = "mysqlnp_empty_password_user";

BOOST_DATA_TEST_CASE_F(any_connection_fixture, regular_password, all_transports)
{
    // Setup
    connect_params params = sample.params;
    params.username = regular_user;
    params.password = regular_passwd;

    // Handshake succeeds
    conn.async_connect(params, as_netresult).validate_no_error();
    check_ssl(conn, sample.expect_ssl);
}

BOOST_DATA_TEST_CASE_F(any_connection_fixture, empty_password, all_transports)
{
    // Setup
    connect_params params = sample.params;
    params.username = empty_user;
    params.password = "";

    // Handshake succeeds
    conn.async_connect(params, as_netresult).validate_no_error();
    check_ssl(conn, sample.expect_ssl);
}

BOOST_DATA_TEST_CASE_F(any_connection_fixture, bad_password, all_transports)
{
    // Setup
    connect_params params = sample.params;
    params.username = regular_user;
    params.password = "bad_password";

    // Handshake fails with the expected error code
    conn.async_connect(params, as_netresult)
        .validate_error_contains(common_server_errc::er_access_denied_error, {"access denied", regular_user});
}

// Spotcheck: mysql_native_password works with old connection
BOOST_FIXTURE_TEST_CASE(tcp_connection_, tcp_connection_fixture)
{
    // Connect succeeds
    conn.async_connect(
            get_tcp_endpoint(),
            connect_params_builder().credentials(regular_user, regular_passwd).build_hparams(),
            as_netresult
    )
        .validate_no_error();
}

BOOST_AUTO_TEST_SUITE_END()  // mysql_native_password

// caching_sha2_password. We acquire a lock on the sha256_mutex
// (dummy table, used as a mutex) to avoid race conditions with other test runs
// (which happens in b2 builds).
// The sha256 cache is shared between all clients.
struct caching_sha2_lock : any_connection_fixture
{
    caching_sha2_lock()
    {
        // Connect
        conn.async_connect(connect_params_builder().credentials("root", "").build(), as_netresult)
            .validate_no_error();

        // Acquire the lock
        results r;
        conn.async_execute("LOCK TABLE sha256_mutex WRITE", r, as_netresult).validate_no_error();

        // The lock is released on fixture destruction, when the connection is closed
    }
};

constexpr const char* regular_user = "csha2p_user";
constexpr const char* regular_passwd = "csha2p_password";
constexpr const char* empty_user = "csha2p_empty_password_user";

BOOST_TEST_DECORATOR(*run_if(&server_features::sha256))
BOOST_AUTO_TEST_SUITE(caching_sha2_password, *boost::unit_test::fixture<caching_sha2_lock>())

static void load_sha256_cache(std::string user, std::string password)
{
    // Connecting as the given user loads the cache
    any_connection_fixture fix;
    fix.conn
        .async_connect(
            connect_params_builder().credentials(std::move(user), std::move(password)).build(),
            as_netresult
        )
        .validate_no_error();
}

static void clear_sha256_cache()
{
    // Issuing a FLUSH PRIVILEGES clears the cache
    any_connection_fixture fix;
    fix.conn.async_connect(connect_params_builder().credentials("root", "").build(), as_netresult)
        .validate_no_error();

    results result;
    fix.conn.async_execute("FLUSH PRIVILEGES", result, as_netresult).validate_no_error();
};

// Cache hit means that we are sending the password hashed, so it is OK to not have SSL for this
BOOST_DATA_TEST_CASE_F(any_connection_fixture, cache_hit, all_transports)
{
    // Setup
    connect_params params = sample.params;
    params.username = regular_user;
    params.password = regular_passwd;
    load_sha256_cache(regular_user, regular_passwd);

    // Handshake succeeds
    conn.async_connect(params, as_netresult).validate_no_error();
    check_ssl(conn, sample.expect_ssl);
}

// Cache miss succeeds only if the underlying transport is secure
BOOST_DATA_TEST_CASE_F(any_connection_fixture, cache_miss_success, secure_transports)
{
    // Setup
    connect_params params = sample.params;
    params.username = regular_user;
    params.password = regular_passwd;
    clear_sha256_cache();

    // Handshake succeeds
    conn.async_connect(params, as_netresult).validate_no_error();
    check_ssl(conn, sample.expect_ssl);
}

// A cache miss would force us send a plaintext password over a non-TLS connection, so we fail
BOOST_FIXTURE_TEST_CASE(cache_miss_error, any_connection_fixture)
{
    // Setup
    connect_params params = connect_params_builder()
                                .ssl(ssl_mode::disable)
                                .credentials(regular_user, regular_passwd)
                                .build();
    clear_sha256_cache();

    // Handshake fails
    conn.async_connect(params, as_netresult).validate_error(client_errc::auth_plugin_requires_ssl);
}

// Empty password users can log in regardless of the SSL usage or cache state
BOOST_DATA_TEST_CASE_F(any_connection_fixture, empty_password_cache_hit, all_transports)
{
    // Setup
    connect_params params = sample.params;
    params.username = empty_user;
    params.password = "";
    load_sha256_cache(empty_user, "");

    // Handshake succeeds
    conn.async_connect(params, as_netresult).validate_no_error();
    check_ssl(conn, sample.expect_ssl);
}

BOOST_DATA_TEST_CASE_F(any_connection_fixture, empty_password_cache_miss, all_transports)
{
    // Setup
    connect_params params = sample.params;
    params.username = empty_user;
    params.password = "";
    clear_sha256_cache();

    // Handshake succeeds
    conn.async_connect(params, as_netresult).validate_no_error();
    check_ssl(conn, sample.expect_ssl);
}

BOOST_FIXTURE_TEST_CASE(bad_password_cache_hit, any_connection_fixture)
{
    // Note: test over non-TLS would return "ssl required"
    auto params = connect_params_builder()
                      .ssl(ssl_mode::require)
                      .credentials(regular_user, "bad_password")
                      .build();
    load_sha256_cache(regular_user, regular_passwd);
    conn.async_connect(params, as_netresult)
        .validate_error_contains(common_server_errc::er_access_denied_error, {"access denied", regular_user});
}

BOOST_FIXTURE_TEST_CASE(bad_password_cache_miss, any_connection_fixture)
{
    // Note: test over non-TLS would return "ssl required"
    auto params = connect_params_builder()
                      .ssl(ssl_mode::require)
                      .credentials(regular_user, "bad_password")
                      .build();
    clear_sha256_cache();
    conn.async_connect(params, as_netresult)
        .validate_error_contains(common_server_errc::er_access_denied_error, {"access denied", regular_user});
}

// Spotcheck: an invalid DB error after cache miss works
BOOST_FIXTURE_TEST_CASE(bad_db_cache_miss, any_connection_fixture)
{
    // Setup
    auto params = connect_params_builder().ssl(ssl_mode::require).database("bad_db").build();
    clear_sha256_cache();

    // Connect fails
    conn.async_connect(params, as_netresult)
        .validate_error(
            common_server_errc::er_dbaccess_denied_error,
            "Access denied for user 'integ_user'@'%' to database 'bad_db'"
        );
}

// Spotcheck: caching_sha2_password works with old connection
BOOST_FIXTURE_TEST_CASE(tcp_ssl_connection_, io_context_fixture)
{
    // Setup
    asio::ssl::context ssl_ctx(asio::ssl::context::tls_client);
    tcp_ssl_connection conn(ctx, ssl_ctx);
    auto params = connect_params_builder().credentials(regular_user, regular_passwd).build_hparams();

    // Connect succeeds
    conn.async_connect(get_tcp_endpoint(), params, as_netresult).validate_no_error();

    // Close succeeds
    conn.async_close(as_netresult).validate_no_error();
}

BOOST_AUTO_TEST_SUITE_END()  // caching_sha2_password

// SSL certificate validation
// This also tests that we can pass a custom ssl::context to connections
BOOST_AUTO_TEST_SUITE(ssl_certificate_validation)

BOOST_AUTO_TEST_CASE(certificate_valid)
{
    // Setup
    asio::ssl::context ssl_ctx(asio::ssl::context::tls_client);
    ssl_ctx.set_verify_mode(boost::asio::ssl::verify_peer);
    ssl_ctx.add_certificate_authority(boost::asio::buffer(CA_PEM));
    any_connection_fixture fix(ssl_ctx);

    // Connect works
    fix.conn.async_connect(connect_params_builder().ssl(ssl_mode::require).build(), as_netresult)
        .validate_no_error();
    check_ssl(fix.conn, true);
}

BOOST_AUTO_TEST_CASE(certificate_invalid)
{
    // Setup
    asio::ssl::context ssl_ctx(asio::ssl::context::tls_client);
    ssl_ctx.set_verify_mode(boost::asio::ssl::verify_peer);
    any_connection_fixture fix(ssl_ctx);

    // Connect fails
    auto err = fix.conn.async_connect(connect_params_builder().ssl(ssl_mode::require).build(), as_netresult)
                   .run()
                   .err;
    BOOST_TEST(err.message().find("certificate verify failed") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(custom_certificate_verification_success)
{
    // Setup
    asio::ssl::context ssl_ctx(asio::ssl::context::tls_client);
    ssl_ctx.set_verify_mode(boost::asio::ssl::verify_peer);
    ssl_ctx.add_certificate_authority(boost::asio::buffer(CA_PEM));
    ssl_ctx.set_verify_callback(boost::asio::ssl::host_name_verification("mysql"));
    any_connection_fixture fix(ssl_ctx);

    // Connect succeeds
    fix.conn.async_connect(connect_params_builder().ssl(ssl_mode::require).build(), as_netresult)
        .validate_no_error();
    check_ssl(fix.conn, true);
}

BOOST_AUTO_TEST_CASE(custom_certificate_verification_error)
{
    // Setup
    asio::ssl::context ssl_ctx(asio::ssl::context::tls_client);
    ssl_ctx.set_verify_mode(boost::asio::ssl::verify_peer);
    ssl_ctx.add_certificate_authority(boost::asio::buffer(CA_PEM));
    ssl_ctx.set_verify_callback(boost::asio::ssl::host_name_verification("host.name"));
    any_connection_fixture fix(ssl_ctx);

    // Connect fails
    auto err = fix.conn.async_connect(connect_params_builder().ssl(ssl_mode::require).build(), as_netresult)
                   .run()
                   .err;
    BOOST_TEST(err.message().find("certificate verify failed") != std::string::npos);
}

// Spotcheck: a custom SSL context can be used with old connections
BOOST_FIXTURE_TEST_CASE(tcp_ssl_connection_, io_context_fixture)
{
    // Setup
    asio::ssl::context ssl_ctx(asio::ssl::context::tls_client);
    ssl_ctx.set_verify_mode(boost::asio::ssl::verify_peer);
    ssl_ctx.add_certificate_authority(boost::asio::buffer(CA_PEM));
    ssl_ctx.set_verify_callback(boost::asio::ssl::host_name_verification("host.name"));
    tcp_ssl_connection conn(ctx, ssl_ctx);
    auto params = connect_params_builder().build_hparams();

    // Connect fails
    auto err = conn.async_connect(get_tcp_endpoint(), params, as_netresult).run().err;
    BOOST_TEST(err.message().find("certificate verify failed") != std::string::npos);
}

BOOST_AUTO_TEST_SUITE_END()  // ssl_certificate_validation

BOOST_AUTO_TEST_SUITE(ssl_mode_)

// All our CI servers support SSL, so enable should behave like required
BOOST_FIXTURE_TEST_CASE(any_enable, any_connection_fixture)
{
    // Setup
    auto params = connect_params_builder().ssl(ssl_mode::enable).build();

    // Connect succeeds
    conn.async_connect(params, as_netresult).validate_no_error();
    check_ssl(conn, true);
}

// connection<>: all ssl modes work as disabled if the stream doesn't support ssl
BOOST_DATA_TEST_CASE_F(
    tcp_connection_fixture,
    non_ssl_stream,
    data::make({ssl_mode::disable, ssl_mode::enable, ssl_mode::require})
)
{
    // Physical connect
    conn.stream().async_connect(get_tcp_endpoint(), as_netresult).validate_no_error_nodiag();

    // Handshake succeeds
    conn.async_handshake(connect_params_builder().ssl(sample).build_hparams(), as_netresult)
        .validate_no_error();
    check_ssl(conn, false);
}

// connection<>: disable can be used to effectively disable SSL
BOOST_AUTO_TEST_CASE(ssl_stream)
{
    struct
    {
        string_view name;
        ssl_mode mode;
        bool expect_ssl;
    } test_cases[] = {
        {"disable", ssl_mode::disable, false},
        {"enable",  ssl_mode::enable,  true },
        {"require", ssl_mode::require, true },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // Setup
            io_context_fixture fix;
            asio::ssl::context ssl_ctx(asio::ssl::context::tls_client);
            tcp_ssl_connection conn(fix.ctx, ssl_ctx);
            auto params = connect_params_builder().ssl(tc.mode).build_hparams();

            // Handshake succeeds
            conn.async_connect(get_tcp_endpoint(), params, as_netresult).validate_no_error();
            check_ssl(conn, tc.expect_ssl);

            // Close succeeds
            conn.async_close(as_netresult).validate_no_error();
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

// Old tcp_ssl_connection, unix_connection, unix_ssl_connection
// can establish and terminate connections, using sync and async fns
BOOST_AUTO_TEST_SUITE(connection_stream_types)

template <class Conn>
struct fixture;

template <>
struct fixture<tcp_ssl_connection> : io_context_fixture
{
    asio::ssl::context ssl_ctx{asio::ssl::context::tls_client};
    tcp_ssl_connection conn{ctx, ssl_ctx};

    using endpoint_type = asio::ip::tcp::endpoint;
    static endpoint_type get_endpoint() { return get_tcp_endpoint(); }
    static bool expect_ssl() { return true; }
};

#ifdef BOOST_ASIO_HAS_LOCAL_SOCKETS
template <>
struct fixture<unix_connection> : io_context_fixture
{
    unix_connection conn{ctx};

    using endpoint_type = asio::local::stream_protocol::endpoint;
    static endpoint_type get_endpoint() { return default_unix_path; }
    static bool expect_ssl() { return false; }
};

template <>
struct fixture<unix_ssl_connection> : io_context_fixture
{
    asio::ssl::context ssl_ctx{asio::ssl::context::tls_client};
    unix_ssl_connection conn{ctx, ssl_ctx};

    using endpoint_type = asio::local::stream_protocol::endpoint;
    static endpoint_type get_endpoint() { return default_unix_path; }
    static bool expect_ssl() { return true; }
};
#endif

template <class Conn>
void do_connect_close_test()
{
    using fixture_type = fixture<Conn>;
    using netmaker_connect = netfun_maker<
        void,
        Conn,
        const typename fixture_type::endpoint_type&,
        const handshake_params&>;
    using netmaker_execute = netfun_maker<void, Conn, const string_view&, results&>;
    using netmaker_close = netfun_maker<void, Conn>;

    struct
    {
        string_view name;
        typename netmaker_connect::signature connect;
        typename netmaker_execute::signature execute;
        typename netmaker_close::signature close;
    } test_cases[] = {
        {"sync",
         netmaker_connect::sync_errc(&Conn::connect),
         netmaker_execute::sync_errc(&Conn::execute),
         netmaker_close::sync_errc(&Conn::close)       },
        {"async",
         netmaker_connect::async_diag(&Conn::async_connect),
         netmaker_execute::async_diag(&Conn::async_execute),
         netmaker_close::async_diag(&Conn::async_close)},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // Setup
            fixture_type fix;

            // Connect
            tc.connect(fix.conn, fix.get_endpoint(), connect_params_builder().build_hparams())
                .validate_no_error();

            // Check whether the connection is using SSL
            check_ssl(fix.conn, fix.expect_ssl());

            // The connection is usable
            results r;
            tc.execute(fix.conn, "SELECT 'abc'", r).validate_no_error();
            BOOST_TEST(r.rows() == makerows(1, "abc"), per_element());

            // Closing succeeds
            tc.close(fix.conn).validate_no_error();
        }
    }
}

template <class Conn>
void do_handshake_quit_test()
{
    using fixture_type = fixture<Conn>;
    using socket_type = typename Conn::stream_type::lowest_layer_type;
    using netmaker_connect = netfun_maker<void, socket_type, const typename fixture_type::endpoint_type&>;
    using netmaker_handshake = netfun_maker<void, Conn, const handshake_params&>;
    using netmaker_execute = netfun_maker<void, Conn, const string_view&, results&>;
    using netmaker_quit = netfun_maker<void, Conn>;

    struct
    {
        string_view name;
        typename netmaker_connect::signature connect;
        typename netmaker_handshake::signature handshake;
        typename netmaker_execute::signature execute;
        typename netmaker_quit::signature quit;
    } test_cases[] = {
        {"sync",
         netmaker_connect::sync_errc_nodiag(&socket_type::connect),
         netmaker_handshake::sync_errc(&Conn::handshake),
         netmaker_execute::sync_errc(&Conn::execute),
         netmaker_quit::sync_errc(&Conn::quit)       },
        {"async",
         netmaker_connect::async_nodiag(&socket_type::async_connect),
         netmaker_handshake::async_diag(&Conn::async_handshake),
         netmaker_execute::async_diag(&Conn::async_execute),
         netmaker_quit::async_diag(&Conn::async_quit)},
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            // Setup
            fixture_type fix;

            // Connect
            tc.connect(fix.conn.stream().lowest_layer(), fix.get_endpoint()).validate_no_error_nodiag();
            tc.handshake(fix.conn, connect_params_builder().build_hparams()).validate_no_error();

            // Check whether the connection uses SSL
            check_ssl(fix.conn, fix.expect_ssl());

            // The connection is usable
            results r;
            tc.execute(fix.conn, "SELECT 'abc'", r).validate_no_error();
            BOOST_TEST(r.rows() == makerows(1, "abc"), per_element());

            // Quitting succeeds
            tc.quit(fix.conn).validate_no_error();
            fix.conn.stream().lowest_layer().close();
        }
    }
}

// tcp_ssl
BOOST_AUTO_TEST_CASE(tcp_ssl_connect_close) { do_connect_close_test<tcp_ssl_connection>(); }
BOOST_AUTO_TEST_CASE(tcp_ssl_handshake_quit) { do_handshake_quit_test<tcp_ssl_connection>(); }

#ifdef BOOST_ASIO_HAS_LOCAL_SOCKETS
// unix
BOOST_TEST_DECORATOR(*run_if(&server_features::unix_sockets))
BOOST_AUTO_TEST_CASE(unix_connection_connect_close) { do_connect_close_test<unix_connection>(); }

BOOST_TEST_DECORATOR(*run_if(&server_features::unix_sockets))
BOOST_AUTO_TEST_CASE(unix_connection_handshake_quit) { do_handshake_quit_test<unix_connection>(); }

// unix ssl
BOOST_TEST_DECORATOR(*run_if(&server_features::unix_sockets))
BOOST_AUTO_TEST_CASE(unix_ssl_connection_connect_close) { do_connect_close_test<unix_ssl_connection>(); }

BOOST_TEST_DECORATOR(*run_if(&server_features::unix_sockets))
BOOST_AUTO_TEST_CASE(unix_ssl_connection_handshake_quit) { do_handshake_quit_test<unix_ssl_connection>(); }
#endif

BOOST_AUTO_TEST_SUITE_END()

// Other handshake tests
BOOST_FIXTURE_TEST_CASE(no_database, any_connection_fixture)
{
    // Connect succeeds
    conn.async_connect(connect_params_builder().database("").build(), as_netresult).validate_no_error();

    // No database selected
    results r;
    conn.async_execute("SELECT DATABASE()", r, as_netresult).validate_no_error();
    BOOST_TEST(r.rows() == makerows(1, nullptr), per_element());
}

BOOST_FIXTURE_TEST_CASE(bad_database, any_connection_fixture)
{
    // Connect fails
    conn.async_connect(connect_params_builder().database("bad_db").build(), as_netresult)
        .validate_error(
            common_server_errc::er_dbaccess_denied_error,
            "Access denied for user 'integ_user'@'%' to database 'bad_db'"
        );
}

BOOST_TEST_DECORATOR(*run_if(&server_features::sha256))
BOOST_FIXTURE_TEST_CASE(unknown_auth_plugin, any_connection_fixture)
{
    // Note: sha256_password is not supported, so it's an unknown plugin to us
    // Setup
    auto params = connect_params_builder()
                      .ssl(ssl_mode::require)
                      .credentials("sha2p_user", "sha2p_password")
                      .build();

    // Connect fails
    conn.async_connect(params, as_netresult).validate_error(client_errc::unknown_auth_plugin);
}

BOOST_FIXTURE_TEST_CASE(bad_user, any_connection_fixture)
{
    // unreliable without SSL. If the default plugin requires SSL
    // (like SHA256), this would fail with 'ssl required'
    // Setup
    auto params = connect_params_builder()
                      .ssl(ssl_mode::require)
                      .credentials("non_existing_user", "bad_password")
                      .build();

    // Connect fails
    auto err = conn.async_connect(params, as_netresult).run().err;
    BOOST_TEST((err.category() == get_common_server_category() || err.category() == get_client_category()));
    BOOST_TEST(err != error_code());  // may be access denied or unknown auth plugin
}

BOOST_AUTO_TEST_SUITE_END()  // test_handshake

}  // namespace
