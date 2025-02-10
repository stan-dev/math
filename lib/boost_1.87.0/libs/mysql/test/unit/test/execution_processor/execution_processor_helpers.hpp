//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_TEST_EXECUTION_PROCESSOR_EXECUTION_PROCESSOR_HELPERS_HPP
#define BOOST_MYSQL_TEST_UNIT_TEST_EXECUTION_PROCESSOR_EXECUTION_PROCESSOR_HELPERS_HPP

#include <boost/mysql/metadata.hpp>
#include <boost/mysql/metadata_collection_view.hpp>

#include <boost/mysql/detail/coldef_view.hpp>

#include <boost/test/unit_test.hpp>

#include "test_common/check_meta.hpp"
#include "test_unit/create_meta.hpp"
#include "test_unit/create_ok.hpp"

// Throughout execution_processor tests we use a set of common values for
// metadata and OK packets. Definitions here to reduce duplication

namespace boost {
namespace mysql {
namespace test {

// Metadata creation
inline detail::coldef_view create_meta_r1_0()
{
    return meta_builder().type(column_type::tinyint).name("ftiny").nullable(false).build_coldef();
}
inline detail::coldef_view create_meta_r1_1()
{
    return meta_builder().type(column_type::varchar).name("fvarchar").nullable(false).build_coldef();
}
inline std::vector<detail::coldef_view> create_meta_r1()
{
    return {
        create_meta_r1_0(),
        create_meta_r1_1(),
    };
}

inline detail::coldef_view create_meta_r2_0()
{
    return meta_builder().type(column_type::bigint).name("fbigint").nullable(false).build_coldef();
}
inline std::vector<detail::coldef_view> create_meta_r2()
{
    return {
        create_meta_r2_0(),
    };
}

inline detail::coldef_view create_meta_r3_0()
{
    return meta_builder().type(column_type::float_).name("ffloat").nullable(false).build_coldef();
}
inline detail::coldef_view create_meta_r3_1()
{
    return meta_builder().type(column_type::double_).name("fdouble").nullable(false).build_coldef();
}
inline detail::coldef_view create_meta_r3_2()
{
    return meta_builder().type(column_type::tinyint).name("ftiny").nullable(false).build_coldef();
}
inline std::vector<detail::coldef_view> create_meta_r3()
{
    return {
        create_meta_r3_0(),
        create_meta_r3_1(),
        create_meta_r3_2(),
    };
}

// Metadata checking
inline void check_meta_r1(metadata_collection_view meta)
{
    check_meta(meta, {column_type::tinyint, column_type::varchar});
}

inline void check_meta_r2(metadata_collection_view meta) { check_meta(meta, {column_type::bigint}); }

inline void check_meta_r3(metadata_collection_view meta)
{
    check_meta(meta, {column_type::float_, column_type::double_, column_type::tinyint});
}

inline void check_meta_empty(metadata_collection_view meta) { BOOST_TEST(meta.size() == 0u); }

// OK packet creation
inline detail::ok_view create_ok_r1(bool more_results = false)
{
    return ok_builder()
        .affected_rows(1)
        .last_insert_id(2)
        .warnings(4)
        .info("Information")
        .more_results(more_results)
        .build();
}

inline detail::ok_view create_ok_r2(bool more_results = false)
{
    return ok_builder()
        .affected_rows(5)
        .last_insert_id(6)
        .warnings(8)
        .info("more_info")
        .more_results(more_results)
        .out_params(true)
        .build();
}

inline detail::ok_view create_ok_r3()
{
    return ok_builder().affected_rows(10).last_insert_id(11).warnings(12).info("").build();
}

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
