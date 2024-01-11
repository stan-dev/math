//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/column_type.hpp>
#include <boost/mysql/mysql_collations.hpp>

#include <boost/mysql/impl/internal/protocol/protocol_field_type.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/printing.hpp"
#include "test_unit/create_meta.hpp"

using namespace boost::mysql::detail;
using namespace boost::mysql::test;
namespace collations = boost::mysql::mysql_collations;
using boost::mysql::column_type;

BOOST_AUTO_TEST_SUITE(test_protocol_field_type)

// Tests edge cases not covered by database_types, where the DB sends
// a protocol_field_type that is supposed not to be sent. Introduced due
// to a bug with recent MariaDB versions that were sending medium_blob only
// if you SELECT'ed TEXT variables
BOOST_AUTO_TEST_CASE(compute_column_type_legacy_types)
{
    struct
    {
        const char* name;
        protocol_field_type proto_type;
        std::uint16_t flags;
        std::uint16_t collation;
        column_type expected;
    } test_cases[] = {
        {"tiny_text",      protocol_field_type::tiny_blob,   0, collations::utf8mb4_general_ci, column_type::text     },
        {"tiny_blob",      protocol_field_type::tiny_blob,   0, collations::binary,             column_type::blob     },
        {"medium_text",
         protocol_field_type::medium_blob,
         0,                                                     collations::utf8mb4_general_ci,
         column_type::text                                                                                            },
        {"medium_blob",    protocol_field_type::medium_blob, 0, collations::binary,             column_type::blob     },
        {"long_text",      protocol_field_type::long_blob,   0, collations::utf8mb4_general_ci, column_type::text     },
        {"long_blob",      protocol_field_type::long_blob,   0, collations::binary,             column_type::blob     },
        {"varchar_string",
         protocol_field_type::varchar,
         0,                                                     collations::utf8mb4_general_ci,
         column_type::varchar                                                                                         },
        {"varchar_binary", protocol_field_type::varchar,     0, collations::binary,             column_type::varbinary},
        {"enum",           protocol_field_type::enum_,       0, collations::utf8mb4_general_ci, column_type::enum_    },
        {"set",            protocol_field_type::set,         0, collations::utf8mb4_general_ci, column_type::set      },
        {"null",           protocol_field_type::null,        0, collations::binary,             column_type::unknown  },
    };

    for (const auto& tc : test_cases)
    {
        BOOST_TEST_CONTEXT(tc.name)
        {
            auto res = compute_column_type(tc.proto_type, tc.flags, tc.collation);
            BOOST_TEST(res == tc.expected);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
