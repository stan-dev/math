//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_MOCK_EXECUTION_PROCESSOR_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_MOCK_EXECUTION_PROCESSOR_HPP

#include <boost/mysql/diagnostics.hpp>
#include <boost/mysql/error_code.hpp>
#include <boost/mysql/field_view.hpp>
#include <boost/mysql/metadata.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/mysql/detail/access.hpp>
#include <boost/mysql/detail/execution_processor/execution_processor.hpp>

#include <boost/config.hpp>
#include <boost/test/unit_test.hpp>

#include <cstddef>

#include "test_unit/fail_count.hpp"

namespace boost {
namespace mysql {
namespace test {

class mock_execution_processor : public detail::execution_processor
{
    struct num_calls_t
    {
        std::size_t reset{};
        std::size_t on_num_meta{};
        std::size_t on_meta{};
        std::size_t on_head_ok_packet{};
        std::size_t on_row_batch_start{};
        std::size_t on_row_batch_finish{};
        std::size_t on_row{};
        std::size_t on_row_ok_packet{};
    };

public:
    // Validate the number of calls safely
    class num_calls_validator
    {
        num_calls_t expected_{};
        num_calls_t actual_{};

    public:
        num_calls_validator(num_calls_t actual) noexcept : actual_(actual) {}

        BOOST_ATTRIBUTE_NODISCARD
        num_calls_validator& reset(std::size_t v) noexcept
        {
            expected_.reset = v;
            return *this;
        }

        BOOST_ATTRIBUTE_NODISCARD
        num_calls_validator& on_num_meta(std::size_t v) noexcept
        {
            expected_.on_num_meta = v;
            return *this;
        }

        BOOST_ATTRIBUTE_NODISCARD
        num_calls_validator& on_meta(std::size_t v) noexcept
        {
            expected_.on_meta = v;
            return *this;
        }

        BOOST_ATTRIBUTE_NODISCARD
        num_calls_validator& on_head_ok_packet(std::size_t v) noexcept
        {
            expected_.on_head_ok_packet = v;
            return *this;
        }

        BOOST_ATTRIBUTE_NODISCARD
        num_calls_validator& on_row_batch_start(std::size_t v) noexcept
        {
            expected_.on_row_batch_start = v;
            return *this;
        }

        BOOST_ATTRIBUTE_NODISCARD
        num_calls_validator& on_row_batch_finish(std::size_t v) noexcept
        {
            expected_.on_row_batch_finish = v;
            return *this;
        }

        BOOST_ATTRIBUTE_NODISCARD
        num_calls_validator& on_row(std::size_t v) noexcept
        {
            expected_.on_row = v;
            return *this;
        }

        BOOST_ATTRIBUTE_NODISCARD
        num_calls_validator& on_row_ok_packet(std::size_t v) noexcept
        {
            expected_.on_row_ok_packet = v;
            return *this;
        }

        void validate()
        {
            BOOST_TEST(expected_.reset == actual_.reset);
            BOOST_TEST(expected_.on_num_meta == actual_.on_num_meta);
            BOOST_TEST(expected_.on_meta == actual_.on_meta);
            BOOST_TEST(expected_.on_head_ok_packet == actual_.on_head_ok_packet);
            BOOST_TEST(expected_.on_row_batch_start == actual_.on_row_batch_start);
            BOOST_TEST(expected_.on_row_batch_finish == actual_.on_row_batch_finish);
            BOOST_TEST(expected_.on_row == actual_.on_row);
            BOOST_TEST(expected_.on_row_ok_packet == actual_.on_row_ok_packet);
        }
    };

    mock_execution_processor() = default;
    std::uint64_t affected_rows() const noexcept { return ok_packet_.affected_rows; }
    std::uint64_t last_insert_id() const noexcept { return ok_packet_.last_insert_id; }
    string_view info() const noexcept { return ok_packet_.info; }
    std::size_t num_meta() const noexcept { return num_meta_; }
    const std::vector<metadata>& meta() const noexcept { return meta_; }
    const std::vector<detail::output_ref>& refs() const noexcept { return refs_; }

    BOOST_ATTRIBUTE_NODISCARD
    num_calls_validator num_calls() noexcept { return num_calls_validator(num_calls_); }

    void set_fail_count(fail_count fc, diagnostics diag = diagnostics())
    {
        fc_ = fc;
        diag_ = std::move(diag);
    }

private:
    // Data
    num_calls_t num_calls_;
    struct
    {
        std::uint64_t affected_rows{};
        std::uint64_t last_insert_id{};
        std::string info;
    } ok_packet_{};
    std::size_t num_meta_{};
    std::vector<metadata> meta_;
    std::vector<detail::output_ref> refs_;
    fail_count fc_;
    diagnostics diag_;

    // Helpers
    error_code maybe_fail(diagnostics& diag)
    {
        auto err = fc_.maybe_fail();
        if (err)
        {
            diag = diag_;
        }
        return err;
    }

    void handle_ok(const detail::ok_view& pack)
    {
        ok_packet_.affected_rows = pack.affected_rows;
        ok_packet_.last_insert_id = pack.last_insert_id;
        ok_packet_.info = pack.info;
    }

protected:
    void reset_impl() noexcept override { ++num_calls_.reset; }
    error_code on_head_ok_packet_impl(const detail::ok_view& pack, diagnostics& diag) override
    {
        ++num_calls_.on_head_ok_packet;
        handle_ok(pack);
        return maybe_fail(diag);
    }
    void on_num_meta_impl(std::size_t num_meta) override
    {
        ++num_calls_.on_num_meta;
        num_meta_ = num_meta;
    }
    error_code on_meta_impl(const detail::coldef_view& coldef, bool is_last, diagnostics& diag) override
    {
        ++num_calls_.on_meta;
        meta_.push_back(detail::access::construct<metadata>(coldef, true));
        return is_last ? maybe_fail(diag) : error_code();
    }
    void on_row_batch_start_impl() override { ++num_calls_.on_row_batch_start; }
    void on_row_batch_finish_impl() override { ++num_calls_.on_row_batch_finish; }
    error_code on_row_impl(span<const std::uint8_t>, const detail::output_ref& ref, std::vector<field_view>&)
        override
    {
        ++num_calls_.on_row;
        refs_.push_back(ref);
        return fc_.maybe_fail();
    }
    error_code on_row_ok_packet_impl(const detail::ok_view& pack) override
    {
        ++num_calls_.on_row_ok_packet;
        handle_ok(pack);
        return fc_.maybe_fail();
    }
};

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
