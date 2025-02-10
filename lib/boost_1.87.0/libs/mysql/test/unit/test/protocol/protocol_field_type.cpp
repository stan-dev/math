//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/column_type.hpp>
#include <boost/mysql/mysql_collations.hpp>

#include <boost/mysql/impl/internal/protocol/impl/protocol_field_type.hpp>

#include <boost/test/unit_test.hpp>

using namespace boost::mysql::detail;
namespace collations = boost::mysql::mysql_collations;
using boost::mysql::column_type;

namespace boost {
namespace mysql {
namespace detail {

// protocol_field_type
inline std::ostream& operator<<(std::ostream& os, protocol_field_type t)
{
    switch (t)
    {
    case protocol_field_type::decimal: return os << "decimal";
    case protocol_field_type::tiny: return os << "tiny";
    case protocol_field_type::short_: return os << "short_";
    case protocol_field_type::long_: return os << "long_";
    case protocol_field_type::float_: return os << "float_";
    case protocol_field_type::double_: return os << "double_";
    case protocol_field_type::null: return os << "null";
    case protocol_field_type::timestamp: return os << "timestamp";
    case protocol_field_type::longlong: return os << "longlong";
    case protocol_field_type::int24: return os << "int24";
    case protocol_field_type::date: return os << "date";
    case protocol_field_type::time: return os << "time";
    case protocol_field_type::datetime: return os << "datetime";
    case protocol_field_type::year: return os << "year";
    case protocol_field_type::varchar: return os << "varchar";
    case protocol_field_type::bit: return os << "bit";
    case protocol_field_type::newdecimal: return os << "newdecimal";
    case protocol_field_type::enum_: return os << "enum_";
    case protocol_field_type::set: return os << "set";
    case protocol_field_type::tiny_blob: return os << "tiny_blob";
    case protocol_field_type::medium_blob: return os << "medium_blob";
    case protocol_field_type::long_blob: return os << "long_blob";
    case protocol_field_type::blob: return os << "blob";
    case protocol_field_type::var_string: return os << "var_string";
    case protocol_field_type::string: return os << "string";
    case protocol_field_type::geometry: return os << "geometry";
    default: return os << "unknown";
    }
}

}  // namespace detail
}  // namespace mysql
}  // namespace boost

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
        {"medium_text",    protocol_field_type::medium_blob, 0, collations::utf8mb4_general_ci, column_type::text
        },
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
