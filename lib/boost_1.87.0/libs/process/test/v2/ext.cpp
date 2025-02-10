// Copyright (c) 2022 Klemens D. Morgenstern
// Copyright (c) 2022 Samuel Venable
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/process/v2/ext/cmd.hpp>
#include <boost/process/v2/ext/cwd.hpp>
#include <boost/process/v2/ext/env.hpp>
#include <boost/process/v2/ext/exe.hpp>
#include <boost/process/v2/pid.hpp>
#include <boost/process/v2/process.hpp>
#include <boost/process/v2/start_dir.hpp>
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(ext)

BOOST_AUTO_TEST_CASE(test_exe)
{
    using boost::unit_test::framework::master_test_suite;

    namespace bp2 = boost::process::v2;
    auto pth = bp2::ext::exe(bp2::current_pid());
    BOOST_CHECK(!pth.empty());
    BOOST_CHECK_EQUAL(bp2::filesystem::canonical(master_test_suite().argv[0]).string(),
                      bp2::filesystem::canonical(pth).string());
}


BOOST_AUTO_TEST_CASE(test_child_exe)
{
    namespace bp2 = boost::process::v2;
    using boost::unit_test::framework::master_test_suite;
    const auto pth = bp2::filesystem::canonical(master_test_suite().argv[0]);

    boost::asio::io_context ctx;
    bp2::process proc(ctx, pth, {"sleep", "10000"});
    BOOST_CHECK_EQUAL(bp2::ext::exe(proc.handle()), pth);
}

BOOST_AUTO_TEST_CASE(cmd)
{
    using boost::unit_test::framework::master_test_suite;

    namespace bp2 = boost::process::v2;
    auto cmd = bp2::ext::cmd(bp2::current_pid());

    // the test framework drops a bunch of args.
    bp2::basic_cstring_ref<typename bp2::shell::char_type> ref(cmd.argv()[0]);
    BOOST_CHECK_EQUAL(
        bp2::detail::conv_string<char>(
            ref.data(), ref.size()
            ), master_test_suite().argv[0]);

    auto cm = cmd.argv() + (cmd.argc() - master_test_suite().argc);
    for (auto i = 1; i < master_test_suite().argc; i++)
    {
      bp2::basic_cstring_ref<typename bp2::shell::char_type> ref(cm[i]);

      BOOST_CHECK_EQUAL(bp2::detail::conv_string<char>(ref.data(), ref.size()),
                        master_test_suite().argv[i]);
    }

}


BOOST_AUTO_TEST_CASE(cmd_exe)
{
    using boost::unit_test::framework::master_test_suite;
    const auto pth =  master_test_suite().argv[0];

    namespace bp2 = boost::process::v2;

    boost::asio::io_context ctx;
    std::vector<std::string> args = {"sleep", "10000", "moar", "args", "  to test "};
    bp2::process proc(ctx, pth, args);
    auto cm = bp2::ext::cmd(proc.handle());

    bp2::basic_cstring_ref<typename bp2::shell::char_type> ref(cm.argv()[0]);
    BOOST_CHECK_EQUAL(bp2::detail::conv_string<char>(ref.data(), ref.size()), pth);

    BOOST_REQUIRE_EQUAL(cm.argc(), args.size() + 1u);
    for (auto i = 0u; i < args.size(); i++)
    {
      ref = cm.argv()[i + 1];

      BOOST_CHECK_EQUAL(bp2::detail::conv_string<char>(ref.data(), ref.size()), args[i]);
    }

}


BOOST_AUTO_TEST_CASE(test_cwd)
{
    namespace bp2 = boost::process::v2;
    auto pth = bp2::ext::cwd(bp2::current_pid()).string();
    if (pth.back() == '\\')
        pth.pop_back();
    BOOST_CHECK_EQUAL(pth, bp2::filesystem::current_path());
}

BOOST_AUTO_TEST_CASE(test_cwd_exe)
{
    using boost::unit_test::framework::master_test_suite;
    namespace bp2 = boost::process::v2;
    const auto pth = bp2::filesystem::absolute(master_test_suite().argv[0]);

    auto tmp = bp2::filesystem::temp_directory_path();

    boost::asio::io_context ctx;
    bp2::process proc(ctx, pth, {"sleep", "10000"},
                      bp2::process_start_dir{tmp});
    auto tt = bp2::ext::cwd(proc.handle()).string();
    if (tt.back() == '\\')
        tt.pop_back();
    BOOST_CHECK_EQUAL(tt, tmp);
    bp2::error_code ec;
    bp2::filesystem::remove(tmp, ec);
}

BOOST_AUTO_TEST_CASE(test_env)
{
    namespace bp2 = boost::process::v2;
    auto env = bp2::ext::env(bp2::current_pid());

    std::size_t e = 0;

    for (const auto & kp : bp2::environment::current())
    {
        auto itr = std::find_if(env.begin(), env.end(),
                                [&](bp2::environment::key_value_pair_view kp_)
                                {
                                   return kp.key() == kp_.key();
                                });
        if (itr != env.end())
        {
            BOOST_CHECK_EQUAL(kp.value(), (*itr).value());
            e++;
        }

    }
    BOOST_CHECK_GT(e, 0u);
}

BOOST_AUTO_TEST_CASE(test_env_exe)
{
    using boost::unit_test::framework::master_test_suite;
    const auto pth =  master_test_suite().argv[0];
    namespace bp2 = boost::process::v2;

    auto tmp = bp2::filesystem::temp_directory_path();

    boost::asio::io_context ctx;

    std::vector<bp2::environment::key_value_pair> new_env;
    {
        auto cr = bp2::environment::current();
        new_env.assign(cr.begin(), cr.end());
    }

    new_env.push_back("FOO=42");
    new_env.push_back("BAR=FOO");

        bp2::process proc(ctx, pth, {"sleep", "10000"},
                      bp2::process_environment(new_env));

    auto env = bp2::ext::env(proc.handle());
    for (const auto & kp : new_env)
    {
        auto itr = std::find_if(env.begin(), env.end(),
                                [&](bp2::environment::key_value_pair_view kp_)
                                {
                                    return kp.key() == kp_.key();
                                });
        BOOST_REQUIRE(itr != env.end());
        BOOST_CHECK_EQUAL(kp.value(), (*itr).value());
    }
}

BOOST_AUTO_TEST_SUITE_END()
