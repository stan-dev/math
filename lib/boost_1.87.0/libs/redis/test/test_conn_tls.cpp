/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <boost/asio/ssl/host_name_verification.hpp>
#include <boost/system/error_code.hpp>
#include <boost/redis/connection.hpp>
#define BOOST_TEST_MODULE conn-tls
#include <boost/test/included/unit_test.hpp>
#include "common.hpp"

namespace net = boost::asio;

using connection = boost::redis::connection;
using boost::redis::request;
using boost::redis::response;
using boost::redis::config;
using boost::system::error_code;

// CA certificate that signed the test server's certificate.
// This is a self-signed CA created for testing purposes.
// This must match tools/tls/ca.crt contents
static constexpr const char* ca_certificate = R"%(-----BEGIN CERTIFICATE-----
MIIFSzCCAzOgAwIBAgIUNd7VUuGK4+ylzCOrmeckg2+TqX8wDQYJKoZIhvcNAQEL
BQAwNTETMBEGA1UECgwKUmVkaXMgVGVzdDEeMBwGA1UEAwwVQ2VydGlmaWNhdGUg
QXV0aG9yaXR5MB4XDTI0MDMzMTE0MjUyM1oXDTM0MDMyOTE0MjUyM1owNTETMBEG
A1UECgwKUmVkaXMgVGVzdDEeMBwGA1UEAwwVQ2VydGlmaWNhdGUgQXV0aG9yaXR5
MIICIjANBgkqhkiG9w0BAQEFAAOCAg8AMIICCgKCAgEA5AMV5V66wt+MM4+oCzH0
xPi++j23p8AOa0o3dxNd4tm5y++gAdKfoxj7oh32ZuYHA5V+sGNEalN/b3GlKXMm
ThdVPSwqOQduny19wrb126ZeQXCfqwgSZQ+rgzaIYpw8/GRRuLDunmsdaR2eiptp
dbv6g6P/aIF6P9mfuekwCC9KBCV6ftqOEnzulNLVw4JjY0rKB9NZqONKVMfWpNyC
zJLCkGmza7BOpybhloZIxGJz033yCjDvIQr9GUWsA5rU9LdUiL+F1W0pWkIel1qo
Evo0EIl3+EOcSSzETI7NPHgnSzNau39ZShV4UBj2lw0DWeNcobeMBQ8ItmqEU6V0
gCEqfUnt10bGIDdmV3D5FKPgvhFvEjQULnblLeLDQ6XDFf+xbGEVjvTzVkLjvyKm
H2D+SKw2O+eDU/0+xhpAf+QsWlm6pmvKWjXI5wK1rh2yssBK2pmY3LuuZCdGrvXb
KX4j/4S9qMr43Hmyoyz0gE5I5rplqot8TvT9O/JsgQYd9fYSvdB+HbqAlJzpBZFl
xbVBXxl0AlDFwQtNMX5ylEQPvYVDKA1M+DTqRTgQKctTfccwvovY3YMV7m5YoODZ
ya2YSBRfQim6VsC+QPYs7p2dk1larIoMMaTaU02oMY+qT2d/eyhWKBv5W9LuowTQ
bWa3ZhWN8lXriPgJOQnZ6iUCAwEAAaNTMFEwHQYDVR0OBBYEFCpEPlClLrgu1zFN
Fmas5G4ybNRJMB8GA1UdIwQYMBaAFCpEPlClLrgu1zFNFmas5G4ybNRJMA8GA1Ud
EwEB/wQFMAMBAf8wDQYJKoZIhvcNAQELBQADggIBAFLl1NZHp0NT5Av4GKmsJFeI
cJOgcIygjR4SBGDAxyPqVpZk0x1q64gJsfOe1ARyI4olQPqO08FZMeB+VBYuqR3S
fEVQZz2FT5U7IVAEZwWHOcWkrrVpEZC6PZktYJ7Yqju6+ic93inoPrHhGNZ5XA/Y
GSfwriWkyWm2SOk35ChFH67MbPWmve8CRAXRmrOCByXwXF87wdqVYZUvH9xDe6WU
snFWXVHr2NA7Re8ZIGp7yJOwwW+CZagepNCPUDwnI0fWOahtOTzonIjq8bfgTZPx
2e7lBuAr9tVMpoeyUytVOlNJDojZAtKOpfMwhAG8ydhk+78aK07VVbnSYVhv7ctU
kkkldqP/S3lBlWo44oOxenwLc9vDQNh64py7eQTD7Qv+TjqAG0ljHIDbVqlkQsgR
pQsu7keG9O1xASSTLZVZN2/alNewpqE/eFRfPM3mtUiTiIZvSxiQnWQMbKofAZH5
HwhVli4RKWRWPqpof4GFNkB8XwfBE+gdlFuWtyg0oRyV3sJ6Zn7E+lUpbQX4CFx3
97vekaFNBchNYMcP3TZ9LwxTx1xOWZ5HHrHyzASG3uz2rqwAsEmdRbmK03KfEQyQ
YpNY718btZ1D6lLino9VMgzaPhUs79bHC64O4ncl7hRclK9qa3KLQdCG1cbIR7G0
2XVYrfsnPHX0CsPDIy7L
-----END CERTIFICATE-----)%";

static config make_tls_config()
{
   config cfg;
   cfg.use_ssl = true;
   cfg.addr.host = get_server_hostname();
   cfg.addr.port = "6380";
   return cfg;
}

BOOST_AUTO_TEST_CASE(ping_internal_ssl_context)
{
   auto const cfg = make_tls_config();
   std::string const in = "Kabuf";

   request req;
   req.push("PING", in);

   response<std::string> resp;

   net::io_context ioc;
   connection conn{ioc};

   // The custom server uses a certificate signed by a CA
   // that is not trusted by default - skip verification.
   conn.next_layer().set_verify_mode(net::ssl::verify_none);

   conn.async_exec(req, resp, [&](error_code ec, auto) {
      BOOST_TEST(ec == std::error_code());
      conn.cancel();
   });

   conn.async_run(cfg, {}, [](auto) { });

   ioc.run();

   BOOST_CHECK_EQUAL(in, std::get<0>(resp).value());
}

BOOST_AUTO_TEST_CASE(ping_custom_ssl_context)
{
   auto const cfg = make_tls_config();
   std::string const in = "Kabuf";

   request req;
   req.push("PING", in);

   response<std::string> resp;

   net::io_context ioc;
   net::ssl::context ctx{boost::asio::ssl::context::tls_client};

   // Configure the SSL context to trust the CA that signed the server's certificate.
   // The test certificate uses "redis" as its common name, regardless of the actual server's hostname
   ctx.add_certificate_authority(net::const_buffer(ca_certificate, std::strlen(ca_certificate)));
   ctx.set_verify_mode(net::ssl::verify_peer);
   ctx.set_verify_callback(net::ssl::host_name_verification("redis"));

   connection conn{ioc, std::move(ctx)};

   conn.async_exec(req, resp, [&](auto ec, auto) {
      BOOST_TEST(ec == std::error_code());
      conn.cancel();
   });

   conn.async_run(cfg, {}, [](auto) { });

   ioc.run();

   BOOST_CHECK_EQUAL(in, std::get<0>(resp).value());
}
