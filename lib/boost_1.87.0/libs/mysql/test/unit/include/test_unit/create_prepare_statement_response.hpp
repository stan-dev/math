//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CREATE_PREPARE_STATEMENT_RESPONSE_HPP
#define BOOST_MYSQL_TEST_UNIT_INCLUDE_TEST_UNIT_CREATE_PREPARE_STATEMENT_RESPONSE_HPP

#include <cstdint>
#include <vector>

namespace boost {
namespace mysql {
namespace test {

class prepare_stmt_response_builder
{
    std::uint8_t seqnum_{};
    std::uint32_t statement_id_{0};
    std::uint16_t num_columns_{0};
    std::uint16_t num_params_{0};

public:
    prepare_stmt_response_builder() = default;

    prepare_stmt_response_builder& seqnum(std::uint8_t v)
    {
        seqnum_ = v;
        return *this;
    }

    prepare_stmt_response_builder& id(std::uint32_t v)
    {
        statement_id_ = v;
        return *this;
    }

    prepare_stmt_response_builder& num_columns(std::uint16_t v)
    {
        num_columns_ = v;
        return *this;
    }

    prepare_stmt_response_builder& num_params(std::uint16_t v)
    {
        num_params_ = v;
        return *this;
    }

    std::vector<std::uint8_t> build() const;
};

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
