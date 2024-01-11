//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

//[example_async_futures

#include <boost/mysql/error_code.hpp>
#include <boost/mysql/handshake_params.hpp>
#include <boost/mysql/results.hpp>
#include <boost/mysql/row_view.hpp>
#include <boost/mysql/statement.hpp>
#include <boost/mysql/tcp_ssl.hpp>

#include <boost/asio/io_context.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/asio/ssl/context.hpp>
#include <boost/asio/use_future.hpp>
#include <boost/system/system_error.hpp>

#include <iostream>
#include <thread>

using boost::asio::use_future;
using boost::mysql::error_code;

void print_employee(boost::mysql::row_view employee)
{
    std::cout << "Employee '" << employee.at(0) << " "   // first_name (string)
              << employee.at(1) << "' earns "            // last_name  (string)
              << employee.at(2) << " dollars yearly\n";  // salary     (double)
}

/**
 * A boost::asio::io_context plus a thread that calls context.run().
 * We encapsulate this here to ensure correct shutdown even in case of
 * error (exception), when we should first reset the work guard, then
 * stop the io_context, and then join the thread. Failing to do so
 * may cause your application to not stop (if the work guard is not
 * reset) or to terminate badly (if the thread is not joined).
 */
class application
{
    boost::asio::io_context ctx_;
    boost::asio::executor_work_guard<boost::asio::io_context::executor_type> guard_;
    std::thread runner_;

public:
    application() : guard_(ctx_.get_executor()), runner_([this] { ctx_.run(); }) {}
    application(const application&) = delete;
    application(application&&) = delete;
    application& operator=(const application&) = delete;
    application& operator=(application&&) = delete;
    ~application()
    {
        guard_.reset();
        runner_.join();
    }
    boost::asio::io_context& context() { return ctx_; }
};

void main_impl(int argc, char** argv)
{
    if (argc != 4 && argc != 5)
    {
        std::cerr << "Usage: " << argv[0] << " <username> <password> <server-hostname> [company-id]\n";
        exit(1);
    }

    // The company_id whose employees we will be listing. This
    // is user-supplied input, and should be treated as untrusted.
    const char* company_id = argc == 5 ? argv[4] : "HGS";

    // Context and connections
    application app;  // boost::asio::io_context and a thread that calls run()
    boost::asio::ssl::context ssl_ctx(boost::asio::ssl::context::tls_client);
    boost::mysql::tcp_ssl_connection conn(app.context(), ssl_ctx);

    // Resolver for hostname resolution
    boost::asio::ip::tcp::resolver resolver(app.context().get_executor());

    // Connection params
    boost::mysql::handshake_params params(
        argv[1],                // username
        argv[2],                // password
        "boost_mysql_examples"  // database to use; leave empty or omit for no database
    );

    /**
     * Hostname resolution.
     * Calling async_resolve triggers the
     * operation, and calling future::get() blocks the current thread until
     * it completes. get() will throw an exception if the operation fails.
     */
    auto endpoints_fut = resolver.async_resolve(
        argv[3],
        boost::mysql::default_port_string,
        boost::asio::use_future
    );
    auto endpoints = endpoints_fut.get();

    // Perform the TCP connect and MySQL handshake.
    std::future<void> fut = conn.async_connect(*endpoints.begin(), params, use_future);
    fut.get();

    // We will be using company_id, which is untrusted user input, so we will use a prepared
    // statement.
    std::future<boost::mysql::statement> stmt_fut = conn.async_prepare_statement(
        "SELECT first_name, last_name, salary FROM employee WHERE company_id = ?",
        use_future
    );
    boost::mysql::statement stmt = stmt_fut.get();

    // Execute the statement
    boost::mysql::results result;
    fut = conn.async_execute(stmt.bind(company_id), result, use_future);
    fut.get();

    // Print employees
    for (boost::mysql::row_view employee : result.rows())
    {
        print_employee(employee);
    }

    // Notify the MySQL server we want to quit, then close the underlying connection.
    conn.async_close(use_future).get();

    // application dtor. stops io_context and then joins the thread
}

int main(int argc, char** argv)
{
    try
    {
        main_impl(argc, argv);
    }
    catch (const std::exception& err)
    {
        std::cerr << "Error: " << err.what() << std::endl;
        return 1;
    }
}

//]
