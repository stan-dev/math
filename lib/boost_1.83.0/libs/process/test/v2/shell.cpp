// Copyright (c) 2022 Klemens D. Morgenstern
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Disable autolinking for unit tests.
#if !defined(BOOST_ALL_NO_LIB)
#define BOOST_ALL_NO_LIB 1
#endif // !defined(BOOST_ALL_NO_LIB)

// Test that header file is self-contained.
#include <boost/process/v2/shell.hpp>
#include <boost/process/v2/process.hpp>
#include <boost/asio/io_context.hpp>

#include <boost/test/unit_test.hpp>

#if defined(BOOST_PROCESS_V2_WINDOWS)
    #define STR(Value) L##Value
    #define STR_VIEW(Value) boost::process::v2::wcstring_ref(STR(Value))
#else
    #define STR(Value) Value
    #define STR_VIEW(Value) boost::process::v2::cstring_ref(STR(Value))
#endif


BOOST_AUTO_TEST_CASE(test_shell_parser)
{
    using boost::process::v2::shell;
    namespace bpv = boost::process::v2;
#if defined(BOOST_PROCESS_V2_POSIX)
    BOOST_CHECK_THROW(shell("foo \""), bpv::system_error);
#endif

    auto sh = shell(STR("foo bar \"foo bar\""));   
    BOOST_CHECK(sh.argc() == 3u);
    BOOST_CHECK(sh.argv()[0] == STR_VIEW("foo"));
    BOOST_CHECK(sh.argv()[1] == STR_VIEW("bar"));
    BOOST_CHECK(sh.argv()[2] == STR_VIEW("foo bar"));

#if defined(BOOST_PROCESS_V2_POSIX)
    auto raw_shell = "sh -c true";
#else
    auto raw_shell = "cmd /c exit 0";
#endif
    sh = shell(raw_shell);

    auto exe = sh.exe();
    BOOST_CHECK(bpv::filesystem::exists(exe));

    boost::asio::io_context ctx;
    bpv::process proc{ctx, exe, sh.args()};

    proc.wait();
    BOOST_CHECK_EQUAL(proc.exit_code(), 0);
}