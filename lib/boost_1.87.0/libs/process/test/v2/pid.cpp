// Copyright (c) 2022 Klemens D. Morgenstern
// Copyright (c) 2022 Samuel Venable
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/process/v2/pid.hpp>
#include <boost/process/v2/process.hpp>

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <thread>
#include <vector>

BOOST_AUTO_TEST_CASE(test_pid)
{
    namespace bp2 = boost::process::v2;
    BOOST_CHECK_NE(bp2::current_pid(), static_cast<bp2::pid_type>(0));

    auto all = bp2::all_pids();
    auto itr = std::find(all.begin(), all.end(), bp2::current_pid());

    BOOST_CHECK_GT(all.size(), 0u);
    BOOST_CHECK(itr != all.end());

}

BOOST_AUTO_TEST_CASE(child_pid)
{
    namespace bp2 = boost::process::v2;
    using boost::unit_test::framework::master_test_suite;
    const auto pth = bp2::filesystem::absolute(master_test_suite().argv[1]);
    std::this_thread::sleep_for(std::chrono::milliseconds(100));

    auto cs = bp2::child_pids(bp2::current_pid());
    boost::asio::io_context ctx;
    bp2::process proc(ctx, pth, {"loop"});
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
    auto c2 = bp2::child_pids(bp2::current_pid());
    BOOST_CHECK_LE(cs.size(), c2.size());
    BOOST_CHECK(std::find(cs.begin(), cs.end(), proc.id()) == cs.end());
    if (!c2.empty())
      BOOST_CHECK(std::find(c2.begin(), c2.end(), proc.id()) != c2.end());
    boost::system::error_code ec;
    proc.terminate(ec);
    if (ec)
      BOOST_CHECK(ec == boost::system::errc::permission_denied);
    else
      proc.wait(ec);

    auto c3 = bp2::child_pids(bp2::current_pid());
    BOOST_CHECK(std::find(c3.begin(), c3.end(), proc.id()) == c3.end());
    BOOST_CHECK_LE(c3.size(), c2.size());
}