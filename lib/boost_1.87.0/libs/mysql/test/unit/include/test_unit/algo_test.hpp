//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_ALGO_TEST_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_ALGO_TEST_HPP

#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>

#include <boost/mysql/detail/next_action.hpp>

#include <boost/mysql/impl/internal/sansio/connection_state_data.hpp>

#include <boost/config.hpp>
#include <boost/core/span.hpp>

#include <cstddef>
#include <cstdint>
#include <vector>

#include "test_common/create_diagnostics.hpp"
#include "test_common/source_location.hpp"

namespace boost {
namespace mysql {
namespace test {

// A type-erased reference to an algorithm
class any_algo_ref
{
    template <class Algo>
    static detail::next_action do_resume(void* self, detail::connection_state_data& st, error_code ec)
    {
        return static_cast<Algo*>(self)->resume(st, ec);
    }

    using fn_t = detail::next_action (*)(void*, detail::connection_state_data&, error_code);

    void* algo_{};
    fn_t fn_{};

public:
    template <class Algo, class = typename std::enable_if<!std::is_same<Algo, any_algo_ref>::value>::type>
    any_algo_ref(Algo& algo) noexcept : algo_(&algo), fn_(&do_resume<Algo>)
    {
    }

    detail::next_action resume(detail::connection_state_data& st, error_code ec)
    {
        return fn_(algo_, st, ec);
    }
};

class BOOST_ATTRIBUTE_NODISCARD algo_test
{
    struct step_t
    {
        detail::next_action_type type;
        std::vector<std::uint8_t> bytes;
        error_code result;
    };

    std::vector<step_t> steps_;

    static void handle_read(detail::connection_state_data& st, const step_t& op);

    detail::next_action run_algo_until_step(
        detail::connection_state_data& st,
        any_algo_ref algo,
        std::size_t num_steps_to_run
    ) const;

    algo_test& add_step(detail::next_action_type act_type, std::vector<std::uint8_t> bytes, error_code ec);

    void check_impl(
        detail::connection_state_data& st,
        any_algo_ref algo,
        const diagnostics& actual_diag,
        error_code expected_ec,
        const diagnostics& expected_diag,
        source_location loc
    ) const;

    std::size_t num_steps() const { return steps_.size(); }

    void check_network_errors_impl(
        detail::connection_state_data& st,
        any_algo_ref algo,
        std::size_t step_number,
        const diagnostics& actual_diag,
        source_location loc
    ) const;

public:
    algo_test() = default;

    BOOST_ATTRIBUTE_NODISCARD
    algo_test& expect_write(std::vector<std::uint8_t> bytes, error_code result = {})
    {
        return add_step(detail::next_action_type::write, std::move(bytes), result);
    }

    BOOST_ATTRIBUTE_NODISCARD
    algo_test& expect_read(std::vector<std::uint8_t> result_bytes)
    {
        return add_step(detail::next_action_type::read, std::move(result_bytes), error_code());
    }

    BOOST_ATTRIBUTE_NODISCARD
    algo_test& expect_read(error_code result) { return add_step(detail::next_action_type::read, {}, result); }

    BOOST_ATTRIBUTE_NODISCARD
    algo_test& expect_ssl_handshake(error_code result = {})
    {
        return add_step(detail::next_action_type::ssl_handshake, {}, result);
    }

    BOOST_ATTRIBUTE_NODISCARD
    algo_test& expect_ssl_shutdown(error_code result = {})
    {
        return add_step(detail::next_action_type::ssl_shutdown, {}, result);
    }

    BOOST_ATTRIBUTE_NODISCARD
    algo_test& expect_close(error_code result = {})
    {
        return add_step(detail::next_action_type::close, {}, result);
    }

    template <class AlgoFixture>
    void check(
        AlgoFixture& fix,
        error_code expected_ec = {},
        const diagnostics& expected_diag = {},
        source_location loc = BOOST_MYSQL_CURRENT_LOCATION
    ) const
    {
        check_impl(fix.st, fix.algo, fix.diag, expected_ec, expected_diag, loc);
    }

    template <class AlgoFixture>
    void check_network_errors(source_location loc = BOOST_MYSQL_CURRENT_LOCATION) const
    {
        for (std::size_t i = 0; i < num_steps(); ++i)
        {
            AlgoFixture fix;
            check_network_errors_impl(fix.st, fix.algo, i, fix.diag, loc);
        }
    }
};

struct algo_fixture_base
{
    static constexpr std::size_t default_max_buffsize = 1024u;

    detail::connection_state_data st;
    diagnostics diag;

    algo_fixture_base(
        diagnostics initial_diag = create_server_diag("Diagnostics not cleared"),
        std::size_t max_buffer_size = default_max_buffsize
    )
        : st(max_buffer_size, max_buffer_size), diag(std::move(initial_diag))
    {
        st.write_buffer.push_back(0xff);  // Check that we clear the write buffer at each step
    }

    algo_fixture_base(std::size_t max_buffer_size) : algo_fixture_base(diagnostics(), max_buffer_size) {}
};

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
