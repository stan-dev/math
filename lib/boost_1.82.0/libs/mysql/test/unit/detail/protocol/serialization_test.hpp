//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_TEST_UNIT_DETAIL_PROTOCOL_SERIALIZATION_TEST_HPP
#define BOOST_MYSQL_TEST_UNIT_DETAIL_PROTOCOL_SERIALIZATION_TEST_HPP

#include <boost/mysql/detail/protocol/constants.hpp>
#include <boost/mysql/detail/protocol/serialization.hpp>

#include <boost/any.hpp>
#include <boost/type_index.hpp>

#include "test_common.hpp"

/**
 * This header defines the required infrastructure for running
 * **serialization tests**. These are based on inspecting real packets
 * with a network analyzer, like Wireshark, to create a "golden file",
 * and verify that our serialization functions generate the same
 * packets as the official MySQL client and server.
 *
 * Each serialization test is defined by an instance of a
 * serialization_sample, which contains a C++ value
 * and its serialized network representation, as a byte array.
 * The interface any_value uses type erasure to represent any C++ value
 * in a generic way, providing access to the required serialization
 * and reporting functions.
 *
 * For each type that can be serialized or deserialized, we define
 * a const serialization_test_spec variable in one of the header
 * files in serialization_test_samples/. This contains a vector of
 * serialization_sample and a serialization_test_type, which defines
 * which of the three following types of test to run:
 *   - serialize: checks serialize() and get_size()
 *   - deserialize: checks deserialize()
 *   - deserialize_space: checks deserialize() under extra bytes and
 *     not enough space conditions. Some messages can't pass these tests,
 *     as their contents depends on message size (e.g. string_eof).
 *
 * Types may run one or more of the above types.
 *
 * All header files in serialization_test_samples/ are included
 * in serialization_test.cpp and the defined variables are used to
 * run the adequate tests.
 */

namespace boost {
namespace mysql {
namespace test {

using detail::int3;
using detail::int_lenenc;
using detail::string_eof;
using detail::string_fixed;
using detail::string_lenenc;
using detail::string_null;

// Helpers for any_value_impl
template <std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<char, N>& v)
{
    return os << string_view(v.data(), N);
}

inline std::ostream& operator<<(std::ostream& os, std::uint8_t value) { return os << +value; }

// Operator == for structs. Very poor efficiency but does the
// job for the purpose of testing
template <class T>
struct struct_member_comparer_op
{
    const T& lhs;
    const T& rhs;
    bool result{};

    struct_member_comparer_op(const T& lhs, const T& rhs) : lhs(lhs), rhs(rhs) {}

    template <class... Types>
    void operator()(const Types&...)
    {
        std::tuple<Types...> lhs_as_tuple;
        std::tuple<Types...> rhs_as_tuple;
        T::apply(lhs, [&lhs_as_tuple](const Types&... args) {
            lhs_as_tuple = std::tuple<Types...>(args...);
        });
        T::apply(rhs, [&rhs_as_tuple](const Types&... args) {
            rhs_as_tuple = std::tuple<Types...>(args...);
        });
        result = lhs_as_tuple == rhs_as_tuple;
    }
};

template <class T>
typename std::enable_if<detail::is_struct_with_fields<T>(), bool>::type operator==(
    const T& lhs,
    const T& rhs
)
{
    struct_member_comparer_op<T> op(lhs, rhs);
    T::apply(lhs, op);
    return op.result;
}

// Operator << for value_holder
template <class T>
std::ostream& operator<<(std::ostream& os, const detail::value_holder<T>& value)
{
    return os << value.value;
}

template <class T>
typename std::enable_if<std::is_enum<T>::value, std::ostream&>::type operator<<(
    std::ostream& os,
    T value
)
{
    return os << boost::typeindex::type_id<T>().pretty_name() << "("
              << static_cast<typename std::underlying_type<T>::type>(value) << ")";
}

// Operator << for structs
struct struct_print_op
{
    std::ostream& os;

    void impl() {}

    template <class T, class... Tail>
    void impl(const T& head, const Tail&... tail)
    {
        os << "    " << head << ",\n";
        impl(tail...);
    }

    template <class... Types>
    void operator()(const Types&... values)
    {
        impl(values...);
    }
};

template <class T>
typename std::enable_if<detail::is_struct_with_fields<T>(), std::ostream&>::type operator<<(
    std::ostream& os,
    const T& value
)
{
    os << boost::typeindex::type_id<T>().pretty_name() << "(\n";
    T::apply(value, struct_print_op{os});
    os << ")\n";
    return os;
}

class any_value
{
public:
    virtual ~any_value() {}
    virtual void serialize(detail::serialization_context& ctx) const = 0;
    virtual std::size_t get_size(const detail::serialization_context& ctx) const = 0;
    virtual detail::deserialize_errc deserialize(detail::deserialization_context& ctx) = 0;
    virtual std::shared_ptr<any_value> default_construct() const = 0;
    virtual bool equals(const any_value& rhs) const = 0;
    virtual void print(std::ostream& os) const = 0;
    virtual std::string type_name() const = 0;

    bool operator==(const any_value& rhs) const { return equals(rhs); }
};
inline std::ostream& operator<<(std::ostream& os, const any_value& value)
{
    value.print(os);
    return os;
}

template <class T>
class any_value_impl : public any_value
{
    T value_;

public:
    any_value_impl(const T& v) : value_(v){};
    void serialize(detail::serialization_context& ctx) const override
    {
        ::boost::mysql::detail::serialize(ctx, value_);
    }
    std::size_t get_size(const detail::serialization_context& ctx) const override
    {
        return ::boost::mysql::detail::get_size(ctx, value_);
    }
    detail::deserialize_errc deserialize(detail::deserialization_context& ctx) override
    {
        return ::boost::mysql::detail::deserialize(ctx, value_);
    }
    std::shared_ptr<any_value> default_construct() const override
    {
        return std::make_shared<any_value_impl<T>>(T{});
    }
    bool equals(const any_value& rhs) const override
    {
        auto typed_value = dynamic_cast<const any_value_impl<T>*>(&rhs);
        return typed_value && (typed_value->value_ == value_);
    }
    void print(std::ostream& os) const override { os << value_; }
    std::string type_name() const override { return boost::typeindex::type_id<T>().pretty_name(); }
};

struct serialization_sample
{
    std::string name;
    std::shared_ptr<any_value> value;
    std::vector<uint8_t> expected_buffer;
    detail::capabilities caps;
    boost::any additional_storage;

    template <class T>
    serialization_sample(
        std::string&& name,
        const T& v,
        std::vector<uint8_t>&& buff,
        std::uint32_t caps = 0,
        boost::any&& storage = {}
    )
        : name(std::move(name)),
          value(std::make_shared<any_value_impl<T>>(v)),
          expected_buffer(std::move(buff)),
          caps(caps),
          additional_storage(std::move(storage))
    {
    }
};

inline std::ostream& operator<<(std::ostream& os, const serialization_sample& input)
{
    return os << "(type=" << input.value->type_name() << ", name=" << input.name << ")";
}

enum class serialization_test_type
{
    serialization,          // only serialization
    deserialization,        // only deserialization
    deserialization_space,  // only deserialization, plus not enough space / extra space tests
    full_no_space,          // serialization + deserialization, but not space tests
    full                    // everything
};

struct serialization_test_spec
{
    serialization_test_type type;
    std::vector<serialization_sample> samples;
};

}  // namespace test
}  // namespace mysql
}  // namespace boost

#endif
