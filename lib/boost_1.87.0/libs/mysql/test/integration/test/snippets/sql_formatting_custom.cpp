//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

// This contains a different definition of formatter<employee>
// Having it separate avoids problems

#include <boost/mysql/format_sql.hpp>
#include <boost/mysql/results.hpp>

#include <boost/test/unit_test.hpp>

#include "test_integration/snippets/get_any_connection.hpp"

using namespace boost::mysql;

namespace {

struct employee
{
    std::string first_name;
    std::string last_name;
    std::string company_id;
};

}  // namespace

namespace boost {
namespace mysql {

//[sql_formatting_formatter_specialization_specifiers
template <>
struct formatter<employee>
{
    // Should we format our employee as an INSERT or an UPDATE?
    // This flag is set by parse and used by format
    bool format_as_update{false};

    const char* parse(const char* it, const char* end)
    {
        // Parse any format specifiers for this type.
        // We recognize the 'u' specifier, which means that we should
        // format the type for an UPDATE statement, instead of an INSERT
        // If no specifier is found, default to INSERT
        // Note that the range may be empty, so we must check for this condition
        if (it != end && *it == 'u')
        {
            // The 'u' specifier is present, record it
            format_as_update = true;

            // Mark the current character as parsed
            ++it;
        }

        // Return a pointer to the first unparsed character.
        // If we didn't manage to parse the entire character range, an error will be emitted.
        // The library already takes care of this.
        return it;
    }

    void format(const employee& emp, format_context_base& ctx) const
    {
        if (format_as_update)
        {
            // This should be suitable to be placed in an UPDATE ... SET statement
            format_sql_to(
                ctx,
                "first_name={}, last_name={}, company_id={}",
                emp.first_name,
                emp.last_name,
                emp.company_id
            );
        }
        else
        {
            // Format only the values, as in INSERT statements
            format_sql_to(ctx, "{}, {}, {}", emp.first_name, emp.last_name, emp.company_id);
        }
    }
};
//]

}  // namespace mysql
}  // namespace boost

namespace {

BOOST_AUTO_TEST_CASE(section_sql_formatting_custom)
{
    auto& conn = test::get_any_connection();
    results r;

    //[sql_formatting_formatter_use_specifiers
    // We can now use the 'u' specifier with employee
    std::string query = boost::mysql::format_sql(
        conn.format_opts().value(),
        "UPDATE employee SET {:u} WHERE id = {}",
        employee{"John", "Doe", "HGS"},
        42
    );

    BOOST_TEST(
        query == "UPDATE employee SET first_name='John', last_name='Doe', company_id='HGS' WHERE id = 42"
    );
    //<-
    conn.execute(query, r);
    //->

    // If no format specifier is provided, we get the old behavior
    query = boost::mysql::format_sql(
        conn.format_opts().value(),
        "INSERT INTO employee (first_name, last_name, company_id) VALUES ({}), ({})",
        employee{"John", "Doe", "HGS"},
        employee{"Rick", "Johnson", "AWC"}
    );

    BOOST_TEST(
        query ==
        "INSERT INTO employee (first_name, last_name, company_id) VALUES "
        "('John', 'Doe', 'HGS'), ('Rick', 'Johnson', 'AWC')"
    );
    //]
    conn.execute(query, r);
}

}  // namespace