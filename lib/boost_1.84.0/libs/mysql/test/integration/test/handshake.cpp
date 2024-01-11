//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/client_errc.hpp>
#include <boost/mysql/common_server_errc.hpp>
#include <boost/mysql/handshake_params.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/tcp_ssl.hpp>

#include <boost/asio/ssl/host_name_verification.hpp>
#include <boost/asio/ssl/verify_mode.hpp>
#include <boost/test/unit_test.hpp>

#include <random>
#include <sstream>

#include "test_integration/common.hpp"
#include "test_integration/er_network_variant.hpp"
#include "test_integration/get_endpoint.hpp"
#include "test_integration/streams.hpp"
#include "test_integration/tcp_network_fixture.hpp"

using namespace boost::mysql::test;

using boost::mysql::client_errc;
using boost::mysql::common_server_errc;
using boost::mysql::error_code;
using boost::mysql::handshake_params;
using boost::mysql::ssl_mode;
using boost::mysql::string_view;
using boost::mysql::tcp_ssl_connection;

namespace {

auto net_samples_ssl = create_network_samples({
    "tcp_ssl_sync_errc",
    "tcp_ssl_async_callback",
});

auto net_samples_nossl = create_network_samples({
    "tcp_sync_errc",
    "tcp_async_callback",
});

auto net_samples_both = create_network_samples({
    "tcp_ssl_sync_errc",
    "tcp_ssl_async_callback",
    "tcp_sync_exc",
    "tcp_async_coroutines",
});

BOOST_AUTO_TEST_SUITE(test_handshake)

struct handshake_fixture : network_fixture
{
    void do_handshake_ok()
    {
        conn->handshake(params).validate_no_error();
        BOOST_TEST(conn->uses_ssl() == var->supports_ssl());
    }

    void do_handshake_ok_ssl()
    {
        params.set_ssl(ssl_mode::require);
        conn->handshake(params).validate_no_error();
        BOOST_TEST(conn->uses_ssl());
    }

    void do_handshake_ok_nossl()
    {
        params.set_ssl(ssl_mode::disable);
        conn->handshake(params).validate_no_error();
        BOOST_TEST(!conn->uses_ssl());
    }
};

// mysql_native_password
BOOST_AUTO_TEST_SUITE(mysql_native_password)

BOOST_MYSQL_NETWORK_TEST(regular_user, handshake_fixture, net_samples_both)
{
    setup_and_physical_connect(sample.net);
    set_credentials("mysqlnp_user", "mysqlnp_password");
    do_handshake_ok();
}

BOOST_MYSQL_NETWORK_TEST(empty_password, handshake_fixture, net_samples_both)
{
    setup_and_physical_connect(sample.net);
    set_credentials("mysqlnp_empty_password_user", "");
    do_handshake_ok();
}

BOOST_MYSQL_NETWORK_TEST(bad_password, handshake_fixture, net_samples_both)
{
    setup_and_physical_connect(sample.net);
    set_credentials("mysqlnp_user", "bad_password");
    conn->handshake(params).validate_error(
        common_server_errc::er_access_denied_error,
        {"access denied", "mysqlnp_user"}
    );
}

BOOST_AUTO_TEST_SUITE_END()  // mysql_native_password

// caching_sha2_password. We create a unique user here to avoid clashes
// with other integration tests running at the same time (which happens in b2 builds).
// We should probably migrate the offending tests to unit tests.
struct caching_sha2_user_creator : tcp_network_fixture
{
    static std::string gen_id()
    {
        constexpr std::size_t len = 10;
        constexpr const char alphanum[] = "0123456789abcdefghijklmnopqrstuvwxyz";

        // Random number generation
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<std::size_t> distrib(0, sizeof(alphanum) - 2);  // NULL terminator

        std::string res;
        res.reserve(len);
        for (std::size_t i = 0; i < len; ++i)
            res += alphanum[distrib(gen)];
        return res;
    }

    static const std::string& regular_username()
    {
        static std::string res = "csha2p_user_" + gen_id();
        return res;
    }

    static const std::string& empty_password_username()
    {
        static std::string res = "csha2p_emptypuser_" + gen_id();
        return res;
    }

    caching_sha2_user_creator()
    {
        std::stringstream query;
        query << "CREATE USER '" << regular_username()
              << "'@'%' IDENTIFIED WITH 'caching_sha2_password' BY 'csha2p_password';"
                 "GRANT ALL PRIVILEGES ON boost_mysql_integtests.*TO '"
              << regular_username()
              << "'@'%';"
                 "CREATE USER '"
              << empty_password_username()
              << "'@'%' IDENTIFIED WITH 'caching_sha2_password' BY '';"
                 "GRANT ALL PRIVILEGES ON boost_mysql_integtests.*TO '"
              << empty_password_username()
              << "'@'%';"
                 "FLUSH PRIVILEGES";

        boost::mysql::results result;
        params.set_username("root");
        params.set_password("");
        params.set_multi_queries(true);
        connect();
        conn.execute(query.str(), result);
    }

    ~caching_sha2_user_creator()
    {
        std::stringstream query;
        query << "DROP USER '" << regular_username()
              << "';"
                 "DROP USER '"
              << empty_password_username() << "'";

        boost::mysql::results result;
        conn.execute(query.str(), result);
    }
};

BOOST_TEST_DECORATOR(*boost::unit_test::label("skip_mysql5"))
BOOST_TEST_DECORATOR(*boost::unit_test::label("skip_mariadb"))
BOOST_AUTO_TEST_SUITE(caching_sha2_password, *boost::unit_test::fixture<caching_sha2_user_creator>())

struct caching_sha2_fixture : handshake_fixture
{
    void load_sha256_cache(string_view user, string_view password)
    {
        tcp_ssl_connection conn(ctx, ssl_ctx);
        conn.connect(get_endpoint<tcp_socket>(), handshake_params(user, password));
        conn.close();
    }

    void clear_sha256_cache()
    {
        tcp_ssl_connection conn(ctx, ssl_ctx);
        boost::mysql::results result;
        conn.connect(get_endpoint<tcp_socket>(), handshake_params("root", ""));
        conn.execute("FLUSH PRIVILEGES", result);
        conn.close();
    }
};

BOOST_MYSQL_NETWORK_TEST(ssl_on_cache_hit, caching_sha2_fixture, net_samples_ssl)
{
    setup_and_physical_connect(sample.net);
    set_credentials(caching_sha2_user_creator::regular_username(), "csha2p_password");
    load_sha256_cache(caching_sha2_user_creator::regular_username(), "csha2p_password");
    do_handshake_ok_ssl();
}

BOOST_MYSQL_NETWORK_TEST(ssl_off_cache_hit, caching_sha2_fixture, net_samples_both)
{
    // As we are sending password hashed, it is OK to not have SSL for this
    setup_and_physical_connect(sample.net);
    set_credentials(caching_sha2_user_creator::regular_username(), "csha2p_password");
    load_sha256_cache(caching_sha2_user_creator::regular_username(), "csha2p_password");
    do_handshake_ok_nossl();
}

BOOST_MYSQL_NETWORK_TEST(ssl_on_cache_miss, caching_sha2_fixture, net_samples_ssl)
{
    setup_and_physical_connect(sample.net);
    set_credentials(caching_sha2_user_creator::regular_username(), "csha2p_password");
    clear_sha256_cache();
    do_handshake_ok_ssl();
}

BOOST_MYSQL_NETWORK_TEST(ssl_off_cache_miss, caching_sha2_fixture, net_samples_both)
{
    // A cache miss would force us send a plaintext password over
    // a non-TLS connection, so we fail
    setup_and_physical_connect(sample.net);
    set_credentials(caching_sha2_user_creator::regular_username(), "csha2p_password");
    clear_sha256_cache();
    params.set_ssl(ssl_mode::disable);
    conn->handshake(params).validate_error(client_errc::auth_plugin_requires_ssl, {});
}

BOOST_MYSQL_NETWORK_TEST(empty_password_ssl_on_cache_hit, caching_sha2_fixture, net_samples_ssl)
{
    setup_and_physical_connect(sample.net);
    set_credentials(caching_sha2_user_creator::empty_password_username(), "");
    load_sha256_cache(caching_sha2_user_creator::empty_password_username(), "");
    do_handshake_ok_ssl();
}

BOOST_MYSQL_NETWORK_TEST(empty_password_ssl_off_cache_hit, caching_sha2_fixture, net_samples_both)
{
    // Empty passwords are allowed over non-TLS connections
    setup_and_physical_connect(sample.net);
    set_credentials(caching_sha2_user_creator::empty_password_username(), "");
    load_sha256_cache(caching_sha2_user_creator::empty_password_username(), "");
    do_handshake_ok_nossl();
}

BOOST_MYSQL_NETWORK_TEST(empty_password_ssl_on_cache_miss, caching_sha2_fixture, net_samples_ssl)
{
    setup_and_physical_connect(sample.net);
    set_credentials(caching_sha2_user_creator::empty_password_username(), "");
    clear_sha256_cache();
    do_handshake_ok_ssl();
}

BOOST_MYSQL_NETWORK_TEST(empty_password_ssl_off_cache_miss, caching_sha2_fixture, net_samples_both)
{
    // Empty passwords are allowed over non-TLS connections
    setup_and_physical_connect(sample.net);
    set_credentials(caching_sha2_user_creator::empty_password_username(), "");
    clear_sha256_cache();
    do_handshake_ok_nossl();
}

BOOST_MYSQL_NETWORK_TEST(bad_password_ssl_on_cache_hit, caching_sha2_fixture, net_samples_ssl)
{
    // Note: test over non-TLS would return "ssl required"
    setup_and_physical_connect(sample.net);
    set_credentials(caching_sha2_user_creator::regular_username(), "bad_password");
    load_sha256_cache(caching_sha2_user_creator::regular_username(), "csha2p_password");
    conn->handshake(params).validate_error(
        common_server_errc::er_access_denied_error,
        {"access denied", caching_sha2_user_creator::regular_username()}
    );
}

BOOST_MYSQL_NETWORK_TEST(bad_password_ssl_on_cache_miss, caching_sha2_fixture, net_samples_ssl)
{
    // Note: test over non-TLS would return "ssl required"
    setup_and_physical_connect(sample.net);
    set_credentials(caching_sha2_user_creator::regular_username(), "bad_password");
    clear_sha256_cache();
    conn->handshake(params).validate_error(
        common_server_errc::er_access_denied_error,
        {"access denied", caching_sha2_user_creator::regular_username()}
    );
}

BOOST_AUTO_TEST_SUITE_END()  // caching_sha2_password

// SSL certificate validation
BOOST_AUTO_TEST_SUITE(ssl_certificate_validation)

// The CA file that signed the server's certificate
constexpr const char CA_PEM[] = R"%(-----BEGIN CERTIFICATE-----
MIIDZzCCAk+gAwIBAgIUWznm2UoxXw3j7HCcp9PpiayTvFQwDQYJKoZIhvcNAQEL
BQAwQjELMAkGA1UEBhMCQVUxEzARBgNVBAgMClNvbWUtU3RhdGUxDjAMBgNVBAoM
BW15c3FsMQ4wDAYDVQQDDAVteXNxbDAgFw0yMDA0MDQxNDMwMjNaGA8zMDE5MDgw
NjE0MzAyM1owQjELMAkGA1UEBhMCQVUxEzARBgNVBAgMClNvbWUtU3RhdGUxDjAM
BgNVBAoMBW15c3FsMQ4wDAYDVQQDDAVteXNxbDCCASIwDQYJKoZIhvcNAQEBBQAD
ggEPADCCAQoCggEBAN0WYdvsDb+a0TxOGPejcwZT0zvTrf921mmDUlrLN1Z0hJ/S
ydgQCSD7Q+6za4lTFZCXcvs52xvvS2gfC0yXyYLCT/jA4RQRxuF+/+w1gDWEbGk0
KzEpsBuKrEIvEaVdoS78SxInnW/aegshdrRRocp4JQ6KHsZgkLTxSwPfYSUmMUo0
cRO0Q/ak3VK8NP13A6ZFvZjrBxjS3cSw9HqilgADcyj1D4EokvfI1C9LrgwgLlZC
XVkjjBqqoMXGGlnXOEK+pm8bU68HM/QvMBkb1Amo8pioNaaYgqJUCP0Ch0iu1nUU
HtsWt6emXv0jANgIW0oga7xcT4MDGN/M+IRWLTECAwEAAaNTMFEwHQYDVR0OBBYE
FNxhaGwf5ePPhzK7yOAKD3VF6wm2MB8GA1UdIwQYMBaAFNxhaGwf5ePPhzK7yOAK
D3VF6wm2MA8GA1UdEwEB/wQFMAMBAf8wDQYJKoZIhvcNAQELBQADggEBAAoeJCAX
IDCFoAaZoQ1niI6Ac/cds8G8It0UCcFGSg+HrZ0YujJxWIruRCUG60Q2OAbEvn0+
uRpTm+4tV1Wt92WFeuRyqkomozx0g4CyfsxGX/x8mLhKPFK/7K9iTXM4/t+xQC4f
J+iRmPVsMKQ8YsHYiWVhlOMH9XJQiqERCB2kOKJCH6xkaF2k0GbM2sGgbS7Z6lrd
fsFTOIVx0VxLVsZnWX3byE9ghnDR5jn18u30Cpb/R/ShxNUGIHqRa4DkM5la6uZX
W1fpSW11JBSUv4WnOO0C2rlIu7UJWOROqZZ0OsybPRGGwagcyff2qVRuI2XFvAMk
OzBrmpfHEhF6NDU=
-----END CERTIFICATE-----
)%";

BOOST_MYSQL_NETWORK_TEST(certificate_valid, handshake_fixture, net_samples_ssl)
{
    // Context changes need to be before setup
    ssl_ctx.set_verify_mode(boost::asio::ssl::verify_peer);
    ssl_ctx.add_certificate_authority(boost::asio::buffer(CA_PEM));
    setup_and_physical_connect(sample.net);
    do_handshake_ok_ssl();
}

BOOST_MYSQL_NETWORK_TEST(certificate_invalid, handshake_fixture, net_samples_ssl)
{
    ssl_ctx.set_verify_mode(boost::asio::ssl::verify_peer);
    setup_and_physical_connect(sample.net);
    auto result = conn->handshake(params);
    BOOST_TEST(result.err.message().find("certificate verify failed") != std::string::npos);
}

BOOST_MYSQL_NETWORK_TEST(custom_certificate_verification_failed, handshake_fixture, net_samples_ssl)
{
    ssl_ctx.set_verify_mode(boost::asio::ssl::verify_peer);
    ssl_ctx.add_certificate_authority(boost::asio::buffer(CA_PEM));
    ssl_ctx.set_verify_callback(boost::asio::ssl::host_name_verification("host.name"));
    setup_and_physical_connect(sample.net);
    auto result = conn->handshake(params);
    BOOST_TEST(result.err.message().find("certificate verify failed") != std::string::npos);
}

BOOST_MYSQL_NETWORK_TEST(custom_certificate_verification_ok, handshake_fixture, net_samples_ssl)
{
    ssl_ctx.set_verify_mode(boost::asio::ssl::verify_peer);
    ssl_ctx.add_certificate_authority(boost::asio::buffer(CA_PEM));
    ssl_ctx.set_verify_callback(boost::asio::ssl::host_name_verification("mysql"));
    setup_and_physical_connect(sample.net);
    do_handshake_ok_ssl();
}

BOOST_AUTO_TEST_SUITE_END()  // ssl_certificate_validation

// Other handshake tests
BOOST_MYSQL_NETWORK_TEST(no_database, handshake_fixture, net_samples_both)
{
    setup_and_physical_connect(sample.net);
    params.set_database("");
    do_handshake_ok();
}

BOOST_TEST_DECORATOR(*boost::unit_test::label("skip_mysql5"))
BOOST_TEST_DECORATOR(*boost::unit_test::label("skip_mariadb"))
BOOST_MYSQL_NETWORK_TEST(unknown_auth_plugin, handshake_fixture, net_samples_ssl)
{
    // Note: sha256_password is not supported, so it's an unknown plugin to us
    setup_and_physical_connect(sample.net);
    set_credentials("sha2p_user", "sha2p_password");
    conn->handshake(params).validate_error(client_errc::unknown_auth_plugin, {});
}

BOOST_MYSQL_NETWORK_TEST(bad_user, handshake_fixture, net_samples_nossl)
{
    // unreliable without SSL. If the default plugin requires SSL
    // (like SHA256), this would fail with 'ssl required'
    setup_and_physical_connect(sample.net);
    set_credentials("non_existing_user", "bad_password");
    conn->handshake(params).validate_any_error();  // may be access denied or unknown auth plugin
}

BOOST_MYSQL_NETWORK_TEST(ssl_disable, handshake_fixture, net_samples_both)
{
    // Both SSL and non-SSL streams will act as non-SSL streams
    setup_and_physical_connect(sample.net);
    params.set_ssl(ssl_mode::disable);
    conn->handshake(params).validate_no_error();
    BOOST_TEST(!conn->uses_ssl());
}

BOOST_MYSQL_NETWORK_TEST(ssl_enable_nonssl_streams, handshake_fixture, net_samples_nossl)
{
    // Ignored by non-ssl streams
    setup_and_physical_connect(sample.net);
    params.set_ssl(ssl_mode::enable);
    conn->handshake(params).validate_no_error();
    BOOST_TEST(!conn->uses_ssl());
}

BOOST_MYSQL_NETWORK_TEST(ssl_enable_ssl_streams, handshake_fixture, net_samples_ssl)
{
    // In all our CI systems, our servers support SSL, so
    // ssl_mode::enable will do the same as ssl_mode::require.
    // We test for this fact.
    setup_and_physical_connect(sample.net);
    params.set_ssl(ssl_mode::enable);
    conn->handshake(params).validate_no_error();
    BOOST_TEST(conn->uses_ssl());
}

BOOST_MYSQL_NETWORK_TEST(ssl_require_nonssl_streams, handshake_fixture, net_samples_nossl)
{
    // Ignored by non-ssl streams
    setup_and_physical_connect(sample.net);
    params.set_ssl(ssl_mode::require);
    conn->handshake(params).validate_no_error();
    BOOST_TEST(!conn->uses_ssl());
}

BOOST_MYSQL_NETWORK_TEST(ssl_require_ssl_streams, handshake_fixture, net_samples_ssl)
{
    setup_and_physical_connect(sample.net);
    params.set_ssl(ssl_mode::require);
    conn->handshake(params).validate_no_error();
    BOOST_TEST(conn->uses_ssl());
}

BOOST_AUTO_TEST_SUITE_END()  // test_handshake

}  // namespace
