/* Copyright (c) 2018-2022 Marcelo Zimbres Silva (mzimbres@gmail.com)
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE.txt)
 */

#include <boost/redis/request.hpp>
#include <boost/redis/response.hpp>
#include <boost/redis/adapter/adapt.hpp>
#include <boost/redis/resp3/parser.hpp>

#define BOOST_TEST_MODULE low level
#include <boost/test/included/unit_test.hpp>

#include <map>
#include <iostream>
#include <optional>
#include <sstream>

// TODO: Test with empty strings.

namespace std
{
auto operator==(boost::redis::ignore_t, boost::redis::ignore_t) noexcept {return true;}
auto operator!=(boost::redis::ignore_t, boost::redis::ignore_t) noexcept {return false;}
}

namespace redis = boost::redis;
namespace resp3 = boost::redis::resp3;
using boost::system::error_code;
using boost::redis::request;
using boost::redis::response;
using boost::redis::generic_response;
using boost::redis::ignore;
using boost::redis::ignore_t;
using boost::redis::adapter::result;
using boost::redis::resp3::parser;
using boost::redis::resp3::parse;
using boost::redis::consume_one;
using boost::redis::error;

using boost::redis::adapter::adapt2;
using node_type = result<resp3::node>;
using vec_type = result<std::vector<std::string>>;
using op_vec_type = result<std::optional<std::vector<std::string>>>;

// Set
using set_type = result<std::set<std::string>>;
using mset_type = result<std::multiset<std::string>>;
using uset_type = result<std::unordered_set<std::string>>;
using muset_type = result<std::unordered_multiset<std::string>>;

// Array
using tuple_int_2 = result<response<int, int>>;
using array_type = result<std::array<int, 3>>;
using array_type2 = result<std::array<int, 1>>;

// Map
using map_type =    result<std::map<std::string, std::string>>;
using mmap_type =   result<std::multimap<std::string, std::string>>;
using umap_type =   result<std::unordered_map<std::string, std::string>>;
using mumap_type =  result<std::unordered_multimap<std::string, std::string>>;
using op_map_type = result<std::optional<std::map<std::string, std::string>>>;
using tuple8_type = result<response<std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string>>;

// Null
using op_type_01 = result<std::optional<bool>>;
using op_type_02 = result<std::optional<int>>;
using op_type_03 = result<std::optional<std::string>>;
using op_type_04 = result<std::optional<std::vector<std::string>>>;
using op_type_05 = result<std::optional<std::list<std::string>>>;
using op_type_06 = result<std::optional<std::map<std::string, std::string>>>;
using op_type_07 = result<std::optional<std::unordered_map<std::string, std::string>>>;
using op_type_08 = result<std::optional<std::set<std::string>>>;
using op_type_09 = result<std::optional<std::unordered_set<std::string>>>;

//-------------------------------------------------------------------

template <class Result>
struct expect {
   std::string in;
   Result expected;
   error_code ec{};
   resp3::type error_type = resp3::type::invalid;
};

template <class Result>
auto make_expected(std::string in, Result expected, error_code ec = {}, resp3::type error_type = resp3::type::invalid)
{
   return expect<Result>{in, expected, ec, error_type};
}

template <class Result>
void test_sync(expect<Result> e)
{
   parser p;
   Result result;
   auto adapter = adapt2(result);
   error_code ec;
   auto const res = parse(p, e.in, adapter, ec);

   BOOST_TEST(res); // None of these tests need more data.

   if (ec) {
      BOOST_CHECK_EQUAL(ec, e.ec);
      return;
   }

   if (result.has_value()) {
      BOOST_TEST(bool(result == e.expected));
      BOOST_CHECK_EQUAL(e.in.size(), p.get_consumed());
   } else {
      BOOST_CHECK_EQUAL(result.error().data_type, e.error_type);
   }
}

template <class Result>
void test_sync2(expect<Result> e)
{
   parser p;
   Result result;
   auto adapter = adapt2(result);
   error_code ec;
   auto const res = parse(p, e.in, adapter, ec);

   BOOST_TEST(res); // None of these tests need more data.
   BOOST_CHECK_EQUAL(ec, e.ec);
}

auto make_blob()
{
   std::string str(100000, 'a');
   str[1000] = '\r';
   str[1001] = '\n';
   return str;

   return str;
}

auto const blob = make_blob();

auto make_blob_string(std::string const& b)
{
   std::string wire;
   wire += '$';
   wire += std::to_string(b.size());
   wire += "\r\n";
   wire += b;
   wire += "\r\n";

   return wire;
}

result<std::optional<int>> op_int_ok = 11;
result<std::optional<bool>> op_bool_ok = true;

// TODO: Test a streamed string that is not finished with a string of
// size 0 but other command comes in.
generic_response streamed_string_e1
{{ {boost::redis::resp3::type::streamed_string, 0, 1, ""}
, {boost::redis::resp3::type::streamed_string_part, 1, 1, "Hell"}
, {boost::redis::resp3::type::streamed_string_part, 1, 1, "o wor"}
, {boost::redis::resp3::type::streamed_string_part, 1, 1, "d"}
, {boost::redis::resp3::type::streamed_string_part, 1, 1, ""}
}};

generic_response streamed_string_e2
{{{resp3::type::streamed_string, 0UL, 1UL, {}}, {resp3::type::streamed_string_part, 1UL, 1UL, {}} }};

generic_response const push_e1a
  {{ {resp3::type::push,          4UL, 0UL, {}}
   , {resp3::type::simple_string, 1UL, 1UL, "pubsub"}
   , {resp3::type::simple_string, 1UL, 1UL, "message"}
   , {resp3::type::simple_string, 1UL, 1UL, "some-channel"}
   , {resp3::type::simple_string, 1UL, 1UL, "some message"}
  }};

generic_response const push_e1b
   {{{resp3::type::push, 0UL, 0UL, {}}}};

generic_response const set_expected1a
   {{{resp3::type::set,            6UL, 0UL, {}}
   , {resp3::type::simple_string,  1UL, 1UL, {"orange"}}
   , {resp3::type::simple_string,  1UL, 1UL, {"apple"}}
   , {resp3::type::simple_string,  1UL, 1UL, {"one"}}
   , {resp3::type::simple_string,  1UL, 1UL, {"two"}}
   , {resp3::type::simple_string,  1UL, 1UL, {"three"}}
   , {resp3::type::simple_string,  1UL, 1UL, {"orange"}}
  }};

mset_type const set_e1f{{"apple", "one", "orange", "orange", "three", "two"}};
uset_type const set_e1c{{"apple", "one", "orange", "three", "two"}};
muset_type const set_e1g{{"apple", "one", "orange", "orange", "three", "two"}};
vec_type const set_e1d = {{"orange", "apple", "one", "two", "three", "orange"}};
op_vec_type const set_expected_1e = set_e1d;

generic_response const array_e1a
{{ {resp3::type::array,       3UL, 0UL, {}}
   , {resp3::type::blob_string, 1UL, 1UL, {"11"}}
   , {resp3::type::blob_string, 1UL, 1UL, {"22"}}
   , {resp3::type::blob_string, 1UL, 1UL, {"3"}}
 }};

result<std::vector<int>> const array_e1b{{11, 22, 3}};
result<std::vector<std::string>> const array_e1c{{"11", "22", "3"}};
result<std::vector<std::string>> const array_e1d{};
generic_response const array_e1e{{{resp3::type::array, 0UL, 0UL, {}}}};
array_type const array_e1f{{11, 22, 3}};
result<std::list<int>> const array_e1g{{11, 22, 3}};
result<std::deque<int>> const array_e1h{{11, 22, 3}};

generic_response const map_expected_1a
{{ {resp3::type::map,         4UL, 0UL, {}}
, {resp3::type::blob_string, 1UL, 1UL, {"key1"}}
, {resp3::type::blob_string, 1UL, 1UL, {"value1"}}
, {resp3::type::blob_string, 1UL, 1UL, {"key2"}}
, {resp3::type::blob_string, 1UL, 1UL, {"value2"}}
, {resp3::type::blob_string, 1UL, 1UL, {"key3"}}
, {resp3::type::blob_string, 1UL, 1UL, {"value3"}}
, {resp3::type::blob_string, 1UL, 1UL, {"key3"}}
, {resp3::type::blob_string, 1UL, 1UL, {"value3"}}
}};

map_type const map_expected_1b
{{ {"key1", "value1"}
, {"key2", "value2"}
, {"key3", "value3"}
}};

umap_type const map_e1g
{{ {"key1", "value1"}
, {"key2", "value2"}
, {"key3", "value3"}
}};

mmap_type const map_e1k
{{ {"key1", "value1"}
, {"key2", "value2"}
, {"key3", "value3"}
, {"key3", "value3"}
}};

mumap_type const map_e1l
{{ {"key1", "value1"}
, {"key2", "value2"}
, {"key3", "value3"}
, {"key3", "value3"}
}};

result<std::vector<std::string>> const map_expected_1c
{{ "key1", "value1"
, "key2", "value2"
, "key3", "value3"
, "key3", "value3"
}};

op_map_type const map_expected_1d = map_expected_1b;

op_vec_type const map_expected_1e = map_expected_1c;

tuple8_type const map_e1f
{ std::string{"key1"}, std::string{"value1"}
, std::string{"key2"}, std::string{"value2"}
, std::string{"key3"}, std::string{"value3"}
, std::string{"key3"}, std::string{"value3"}
};

generic_response const attr_e1a
{{ {resp3::type::attribute,     1UL, 0UL, {}}
   , {resp3::type::simple_string, 1UL, 1UL, "key-popularity"}
   , {resp3::type::map,           2UL, 1UL, {}}
   , {resp3::type::blob_string,   1UL, 2UL, "a"}
   , {resp3::type::doublean,      1UL, 2UL, "0.1923"}
   , {resp3::type::blob_string,   1UL, 2UL, "b"}
   , {resp3::type::doublean,      1UL, 2UL, "0.0012"}
}  };

generic_response const attr_e1b
   {{{resp3::type::attribute, 0UL, 0UL, {}} }};

#define S01a  "#11\r\n"
#define S01b  "#f\r\n"
#define S01c  "#t\r\n"
#define S01d  "#\r\n"

#define S02a  "$?\r\n;0\r\n"
#define S02b  "$?\r\n;4\r\nHell\r\n;5\r\no wor\r\n;1\r\nd\r\n;0\r\n"
#define S02c  "$?\r\n;b\r\nHell\r\n;5\r\no wor\r\n;1\r\nd\r\n;0\r\n"
#define S02d  "$?\r\n;d\r\nHell\r\n;5\r\no wor\r\n;1\r\nd\r\n;0\r\n"

#define S03a  "%11\r\n"
#define S03b  "%4\r\n$4\r\nkey1\r\n$6\r\nvalue1\r\n$4\r\nkey2\r\n$6\r\nvalue2\r\n$4\r\nkey3\r\n$6\r\nvalue3\r\n$4\r\nkey3\r\n$6\r\nvalue3\r\n"
#define S03c  "%0\r\n"
#define S03d  "%rt\r\n$4\r\nkey1\r\n$6\r\nvalue1\r\n$4\r\nkey2\r\n$6\r\nvalue2\r\n"

#define S04a  "*1\r\n:11\r\n"
#define S04b  "*3\r\n$2\r\n11\r\n$2\r\n22\r\n$1\r\n3\r\n"
#define S04c  "*1\r\n" S03b
#define S04d  "*1\r\n" S09a
#define S04e  "*3\r\n$2\r\n11\r\n$2\r\n22\r\n$1\r\n3\r\n"
#define S04f  "*1\r\n*1\r\n$2\r\nab\r\n"
#define S04g  "*1\r\n*1\r\n*1\r\n*1\r\n*1\r\n*1\r\na\r\n"
#define S04h  "*0\r\n"
#define S04i  "*3\r\n$2\r\n11\r\n$2\r\n22\r\n$1\r\n3\r\n"

#define S05a  ":-3\r\n"
#define S05b  ":11\r\n"
#define s05c  ":3\r\n"
#define S05d  ":adf\r\n"
#define S05e  ":\r\n"

#define S06a  "_\r\n"

#define S07a  ">4\r\n+pubsub\r\n+message\r\n+some-channel\r\n+some message\r\n"
#define S07b  ">0\r\n"

#define S08a  "|1\r\n+key-popularity\r\n%2\r\n$1\r\na\r\n,0.1923\r\n$1\r\nb\r\n,0.0012\r\n"
#define S08b  "|0\r\n"

#define S09a  "~6\r\n+orange\r\n+apple\r\n+one\r\n+two\r\n+three\r\n+orange\r\n"
#define S09b  "~0\r\n"

#define S10a  "-Error\r\n"
#define S10b  "-\r\n"

#define S11a  ",1.23\r\n"
#define S11b  ",inf\r\n"
#define S11c  ",-inf\r\n" 
#define S11d  ",1.23\r\n"
#define S11e  ",er\r\n"
#define S11f  ",\r\n"

#define S12a  "!21\r\nSYNTAX invalid syntax\r\n"
#define S12b  "!0\r\n\r\n"
#define S12c  "!3\r\nfoo\r\n"

#define S13a  "=15\r\ntxt:Some string\r\n"
#define S13b  "=0\r\n\r\n"

#define S14a  "(3492890328409238509324850943850943825024385\r\n"
#define S14b  "(\r\n"

#define S15a  "+OK\r\n"
#define S15b  "+\r\n"

#define S16a  "s11\r\n"

#define S17a  "$l\r\nhh\r\n"
#define S17b  "$2\r\nhh\r\n"
#define S18c  "$26\r\nhhaa\aaaa\raaaaa\r\naaaaaaaaaa\r\n"
#define S18d  "$0\r\n\r\n" 

#define NUMBER_TEST_CONDITIONS(test) \
   test(make_expected(S01a, result<std::optional<bool>>{}, boost::redis::error::unexpected_bool_value)); \
   test(make_expected(S01b, result<bool>{{false}})); \
   test(make_expected(S01b, node_type{{resp3::type::boolean, 1UL, 0UL, {"f"}}})); \
   test(make_expected(S01c, result<bool>{{true}})); \
   test(make_expected(S01c, node_type{{resp3::type::boolean, 1UL, 0UL, {"t"}}})); \
   test(make_expected(S01c, op_bool_ok)); \
   test(make_expected(S01c, result<std::map<int, int>>{}, boost::redis::error::expects_resp3_map)); \
   test(make_expected(S01c, result<std::set<int>>{}, boost::redis::error::expects_resp3_set)); \
   test(make_expected(S01c, result<std::unordered_map<int, int>>{}, boost::redis::error::expects_resp3_map)); \
   test(make_expected(S01c, result<std::unordered_set<int>>{}, boost::redis::error::expects_resp3_set)); \
   test(make_expected(S02a, streamed_string_e2)); \
   test(make_expected(S03a, result<int>{}, boost::redis::error::expects_resp3_simple_type));\
   test(make_expected(S03a, result<std::optional<int>>{}, boost::redis::error::expects_resp3_simple_type));; \
   test(make_expected(S02b, result<int>{}, boost::redis::error::not_a_number)); \
   test(make_expected(S02b, result<std::string>{std::string{"Hello word"}})); \
   test(make_expected(S02b, streamed_string_e1)); \
   test(make_expected(S02c, result<std::string>{}, boost::redis::error::not_a_number)); \
   test(make_expected(S05a, node_type{{resp3::type::number, 1UL, 0UL, {"-3"}}})); \
   test(make_expected(S05b, result<int>{11})); \
   test(make_expected(S05b, op_int_ok)); \
   test(make_expected(S05b, result<std::list<std::string>>{}, boost::redis::error::expects_resp3_aggregate)); \
   test(make_expected(S05b, result<std::map<std::string, std::string>>{}, boost::redis::error::expects_resp3_map)); \
   test(make_expected(S05b, result<std::set<std::string>>{}, boost::redis::error::expects_resp3_set)); \
   test(make_expected(S05b, result<std::unordered_map<std::string, std::string>>{}, boost::redis::error::expects_resp3_map)); \
   test(make_expected(S05b, result<std::unordered_set<std::string>>{}, boost::redis::error::expects_resp3_set)); \
   test(make_expected(s05c, array_type2{}, boost::redis::error::expects_resp3_aggregate));\
   test(make_expected(s05c, node_type{{resp3::type::number, 1UL, 0UL, {"3"}}}));\
   test(make_expected(S06a, op_type_01{})); \
   test(make_expected(S06a, op_type_02{}));\
   test(make_expected(S06a, op_type_03{}));\
   test(make_expected(S06a, op_type_04{}));\
   test(make_expected(S06a, op_type_05{}));\
   test(make_expected(S06a, op_type_06{}));\
   test(make_expected(S06a, op_type_07{}));\
   test(make_expected(S06a, op_type_08{}));\
   test(make_expected(S06a, op_type_09{}));\
   test(make_expected(S07a, push_e1a)); \
   test(make_expected(S07b, push_e1b)); \
   test(make_expected(S04b, map_type{}, boost::redis::error::expects_resp3_map));\
   test(make_expected(S03b, map_e1f));\
   test(make_expected(S03b, map_e1g));\
   test(make_expected(S03b, map_e1k));\
   test(make_expected(S03b, map_expected_1a));\
   test(make_expected(S03b, map_expected_1b));\
   test(make_expected(S03b, map_expected_1c));\
   test(make_expected(S03b, map_expected_1d));\
   test(make_expected(S03b, map_expected_1e));\
   test(make_expected(S08a, attr_e1a)); \
   test(make_expected(S08b, attr_e1b)); \
   test(make_expected(S04e, array_e1a));\
   test(make_expected(S04e, array_e1b));\
   test(make_expected(S04e, array_e1c));\
   test(make_expected(S04e, array_e1f));\
   test(make_expected(S04e, array_e1g));\
   test(make_expected(S04e, array_e1h));\
   test(make_expected(S04e, array_type2{}, boost::redis::error::incompatible_size));\
   test(make_expected(S04e, tuple_int_2{}, boost::redis::error::incompatible_size));\
   test(make_expected(S04f, array_type2{}, boost::redis::error::nested_aggregate_not_supported));\
   test(make_expected(S04g, generic_response{}, boost::redis::error::exceeeds_max_nested_depth));\
   test(make_expected(S04h, array_e1d));\
   test(make_expected(S04h, array_e1e));\
   test(make_expected(S04i, set_type{}, boost::redis::error::expects_resp3_set)); \
   test(make_expected(S09a, set_e1c)); \
   test(make_expected(S09a, set_e1d)); \
   test(make_expected(S09a, set_e1f)); \
   test(make_expected(S09a, set_e1g)); \
   test(make_expected(S09a, set_expected1a)); \
   test(make_expected(S09a, set_expected_1e)); \
   test(make_expected(S09a, set_type{{"apple", "one", "orange", "three", "two"}})); \
   test(make_expected(S09b, generic_response{{{resp3::type::set,  0UL, 0UL, {}}}})); \
   test(make_expected(S03c, map_type{}));\
   test(make_expected(S11a, node_type{{resp3::type::doublean, 1UL, 0UL, {"1.23"}}}));\
   test(make_expected(S11b, node_type{{resp3::type::doublean, 1UL, 0UL, {"inf"}}}));\
   test(make_expected(S11c, node_type{{resp3::type::doublean, 1UL, 0UL, {"-inf"}}}));\
   test(make_expected(S11d, result<double>{{1.23}}));\
   test(make_expected(S11e, result<double>{{0}}, boost::redis::error::not_a_double));\
   test(make_expected(S13a, node_type{{resp3::type::verbatim_string, 1UL, 0UL, {"txt:Some string"}}}));\
   test(make_expected(S13b, node_type{{resp3::type::verbatim_string, 1UL, 0UL, {}}}));\
   test(make_expected(S14a, node_type{{resp3::type::big_number, 1UL, 0UL, {"3492890328409238509324850943850943825024385"}}}));\
   test(make_expected(S14b, result<int>{}, boost::redis::error::empty_field));\
   test(make_expected(S15a, result<std::optional<std::string>>{{"OK"}}));\
   test(make_expected(S15a, result<std::string>{{"OK"}}));\
   test(make_expected(S15b, result<std::optional<std::string>>{""}));\
   test(make_expected(S15b, result<std::string>{{""}}));\
   test(make_expected(S16a, result<int>{}, boost::redis::error::invalid_data_type));\
   test(make_expected(S05d, result<int>{11}, boost::redis::error::not_a_number));\
   test(make_expected(S03d, map_type{}, boost::redis::error::not_a_number));\
   test(make_expected(S02d, result<std::string>{}, boost::redis::error::not_a_number));\
   test(make_expected(S17a, result<std::string>{}, boost::redis::error::not_a_number));\
   test(make_expected(S05e, result<int>{}, boost::redis::error::empty_field));\
   test(make_expected(S01d, result<std::optional<bool>>{}, boost::redis::error::empty_field));\
   test(make_expected(S11f, result<std::string>{}, boost::redis::error::empty_field));\
   test(make_expected(S17b, node_type{{resp3::type::blob_string, 1UL, 0UL, {"hh"}}}));\
   test(make_expected(S18c, node_type{{resp3::type::blob_string, 1UL, 0UL, {"hhaa\aaaa\raaaaa\r\naaaaaaaaaa"}}}));\
   test(make_expected(S18d, node_type{{resp3::type::blob_string, 1UL, 0UL, {}}}));\
   test(make_expected(make_blob_string(blob), node_type{{resp3::type::blob_string, 1UL, 0UL, {blob}}}));\
   test(make_expected(S04a, result<std::vector<int>>{{11}})); \
   test(make_expected(S04d, result<response<std::unordered_set<std::string>>>{response<std::unordered_set<std::string>>{{set_e1c}}})); \
   test(make_expected(S04c, result<response<std::map<std::string, std::string>>>{response<std::map<std::string, std::string>>{{map_expected_1b}}}));\
   test(make_expected(S03b, map_e1l));\
   test(make_expected(S06a, result<int>{0}, {}, resp3::type::null)); \
   test(make_expected(S06a, map_type{}, {}, resp3::type::null));\
   test(make_expected(S06a, array_type{}, {}, resp3::type::null));\
   test(make_expected(S06a, result<std::list<int>>{}, {}, resp3::type::null));\
   test(make_expected(S06a, result<std::vector<int>>{}, {}, resp3::type::null));\
   test(make_expected(S10a, result<ignore_t>{}, boost::redis::error::resp3_simple_error)); \
   test(make_expected(S10a, node_type{{resp3::type::simple_error, 1UL, 0UL, {"Error"}}}, {}, resp3::type::simple_error)); \
   test(make_expected(S10b, node_type{{resp3::type::simple_error, 1UL, 0UL, {""}}}, {}, resp3::type::simple_error)); \
   test(make_expected(S12a, node_type{{resp3::type::blob_error, 1UL, 0UL, {"SYNTAX invalid syntax"}}}, {}, resp3::type::blob_error));\
   test(make_expected(S12b, node_type{{resp3::type::blob_error, 1UL, 0UL, {}}}, {}, resp3::type::blob_error));\
   test(make_expected(S12c, result<ignore_t>{}, boost::redis::error::resp3_blob_error));\

BOOST_AUTO_TEST_CASE(sansio)
{
   NUMBER_TEST_CONDITIONS(test_sync)
}

BOOST_AUTO_TEST_CASE(ignore_adapter_simple_error)
{
   test_sync2(make_expected(S10a, ignore, boost::redis::error::resp3_simple_error));
}

BOOST_AUTO_TEST_CASE(ignore_adapter_blob_error)
{
   test_sync2(make_expected(S12a, ignore, boost::redis::error::resp3_blob_error));
}

BOOST_AUTO_TEST_CASE(ignore_adapter_no_error)
{
   test_sync2(make_expected(S05b, ignore));
}

//-----------------------------------------------------------------------------------
void check_error(char const* name, boost::redis::error ev)
{
   auto const ec = boost::redis::make_error_code(ev);
   auto const& cat = ec.category();
   BOOST_TEST(std::string(ec.category().name()) == name);
   BOOST_TEST(!ec.message().empty());
   BOOST_TEST(cat.equivalent(
      static_cast<std::underlying_type<boost::redis::error>::type>(ev),
          ec.category().default_error_condition(
              static_cast<std::underlying_type<boost::redis::error>::type>(ev))));
   BOOST_TEST(cat.equivalent(ec,
      static_cast<std::underlying_type<boost::redis::error>::type>(ev)));
}

BOOST_AUTO_TEST_CASE(cover_error)
{
   check_error("boost.redis", boost::redis::error::invalid_data_type);
   check_error("boost.redis", boost::redis::error::not_a_number);
   check_error("boost.redis", boost::redis::error::exceeeds_max_nested_depth);
   check_error("boost.redis", boost::redis::error::unexpected_bool_value);
   check_error("boost.redis", boost::redis::error::empty_field);
   check_error("boost.redis", boost::redis::error::expects_resp3_simple_type);
   check_error("boost.redis", boost::redis::error::expects_resp3_aggregate);
   check_error("boost.redis", boost::redis::error::expects_resp3_map);
   check_error("boost.redis", boost::redis::error::expects_resp3_set);
   check_error("boost.redis", boost::redis::error::nested_aggregate_not_supported);
   check_error("boost.redis", boost::redis::error::resp3_simple_error);
   check_error("boost.redis", boost::redis::error::resp3_blob_error);
   check_error("boost.redis", boost::redis::error::incompatible_size);
   check_error("boost.redis", boost::redis::error::not_a_double);
   check_error("boost.redis", boost::redis::error::resp3_null);
   check_error("boost.redis", boost::redis::error::not_connected);
   check_error("boost.redis", boost::redis::error::resolve_timeout);
   check_error("boost.redis", boost::redis::error::connect_timeout);
   check_error("boost.redis", boost::redis::error::pong_timeout);
   check_error("boost.redis", boost::redis::error::ssl_handshake_timeout);
   check_error("boost.redis", boost::redis::error::sync_receive_push_failed);
   check_error("boost.redis", boost::redis::error::incompatible_node_depth);
}

std::string get_type_as_str(boost::redis::resp3::type t)
{
   std::ostringstream ss;
   ss << t;
   return ss.str();
}

BOOST_AUTO_TEST_CASE(type_string)
{
   BOOST_TEST(!get_type_as_str(boost::redis::resp3::type::array).empty());
   BOOST_TEST(!get_type_as_str(boost::redis::resp3::type::push).empty());
   BOOST_TEST(!get_type_as_str(boost::redis::resp3::type::set).empty());
   BOOST_TEST(!get_type_as_str(boost::redis::resp3::type::map).empty());
   BOOST_TEST(!get_type_as_str(boost::redis::resp3::type::attribute).empty());
   BOOST_TEST(!get_type_as_str(boost::redis::resp3::type::simple_string).empty());
   BOOST_TEST(!get_type_as_str(boost::redis::resp3::type::simple_error).empty());
   BOOST_TEST(!get_type_as_str(boost::redis::resp3::type::number).empty());
   BOOST_TEST(!get_type_as_str(boost::redis::resp3::type::doublean).empty());
   BOOST_TEST(!get_type_as_str(boost::redis::resp3::type::boolean).empty());
   BOOST_TEST(!get_type_as_str(boost::redis::resp3::type::big_number).empty());
   BOOST_TEST(!get_type_as_str(boost::redis::resp3::type::null).empty());
   BOOST_TEST(!get_type_as_str(boost::redis::resp3::type::blob_error).empty());
   BOOST_TEST(!get_type_as_str(boost::redis::resp3::type::verbatim_string).empty());
   BOOST_TEST(!get_type_as_str(boost::redis::resp3::type::blob_string).empty());
   BOOST_TEST(!get_type_as_str(boost::redis::resp3::type::streamed_string_part).empty());
   BOOST_TEST(!get_type_as_str(boost::redis::resp3::type::invalid).empty());
}

BOOST_AUTO_TEST_CASE(type_convert)
{
   using boost::redis::resp3::to_code;
   using boost::redis::resp3::to_type;
   using boost::redis::resp3::type;

#define CHECK_CASE(A) BOOST_CHECK_EQUAL(to_type(to_code(type::A)), type::A);
   CHECK_CASE(array);
   CHECK_CASE(push);
   CHECK_CASE(set);
   CHECK_CASE(map);
   CHECK_CASE(attribute);
   CHECK_CASE(simple_string);
   CHECK_CASE(simple_error);
   CHECK_CASE(number);
   CHECK_CASE(doublean);
   CHECK_CASE(boolean);
   CHECK_CASE(big_number);
   CHECK_CASE(null);
   CHECK_CASE(blob_error);
   CHECK_CASE(verbatim_string);
   CHECK_CASE(blob_string);
   CHECK_CASE(streamed_string_part);
#undef CHECK_CASE
}

BOOST_AUTO_TEST_CASE(adapter)
{
   using boost::redis::adapter::boost_redis_adapt;
   using resp3::type;

   error_code ec;

   response<std::string, int, ignore_t> resp;

   auto f = boost_redis_adapt(resp);
   f(0, resp3::basic_node<std::string_view>{type::simple_string, 1, 0, "Hello"}, ec);
   f(1, resp3::basic_node<std::string_view>{type::number, 1, 0, "42"}, ec);

   BOOST_CHECK_EQUAL(std::get<0>(resp).value(), "Hello");
   BOOST_TEST(!ec);

   BOOST_CHECK_EQUAL(std::get<1>(resp).value(), 42);
   BOOST_TEST(!ec);
}

// TODO: This was an experiment, I will resume implementing this
// later.
BOOST_AUTO_TEST_CASE(adapter_as)
{
   result<std::set<std::string>> set;
   auto adapter = adapt2(set);

   for (auto const& e: set_expected1a.value()) {
      error_code ec;
      adapter(e, ec);
   }
}

BOOST_AUTO_TEST_CASE(cancel_one_1)
{
   auto resp = push_e1a;
   BOOST_TEST(resp.has_value());

   consume_one(resp);
   BOOST_TEST(resp.value().empty());
}

BOOST_AUTO_TEST_CASE(cancel_one_empty)
{
   generic_response resp;
   BOOST_TEST(resp.has_value());

   consume_one(resp);
   BOOST_TEST(resp.value().empty());
}

BOOST_AUTO_TEST_CASE(cancel_one_has_error)
{
   generic_response resp = boost::redis::adapter::error{resp3::type::simple_string, {}};
   BOOST_TEST(resp.has_error());

   consume_one(resp);
   BOOST_TEST(resp.has_error());
}

BOOST_AUTO_TEST_CASE(cancel_one_has_does_not_consume_past_the_end)
{
   auto resp = push_e1a;
   BOOST_TEST(resp.has_value());
   resp.value().insert(
      std::cend(resp.value()),
      std::cbegin(push_e1a.value()),
      std::cend(push_e1a.value()));

   consume_one(resp);

   BOOST_CHECK_EQUAL(resp.value().size(), push_e1a.value().size());
}

BOOST_AUTO_TEST_CASE(cancel_one_incompatible_depth)
{
   auto resp = streamed_string_e1;
   BOOST_TEST(resp.has_value());

   error_code ec;
   consume_one(resp, ec);

   error_code expected = error::incompatible_node_depth;
   BOOST_CHECK_EQUAL(ec, expected);

   BOOST_CHECK_EQUAL(resp.value().size(), push_e1a.value().size());
}
