//
// Copyright (c) 2023 Dmitry Arkhipov (grisumbras@yandex.ru)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/json
//

//[example_use_allocator

/*
    This example uses a context that stores an allocator to create sequences during conversions
*/

#include <boost/container/pmr/monotonic_buffer_resource.hpp>
#include <boost/container/pmr/vector.hpp>
#include <boost/json.hpp>
#include <boost/system/errc.hpp>

using namespace boost::json;
using namespace boost::container;

template< class Alloc >
struct use_allocator_t
{
    Alloc allocator;
};

template< class Alloc >
use_allocator_t< Alloc >
use_allocator( Alloc alloc ) noexcept
{
    return { alloc };
}

template<
    class T,
    class Alloc,
    class FullContext,
    class = typename std::enable_if<
        is_sequence_like<T>::value &&
        std::uses_allocator<T, Alloc>::value
        >::type >
boost::system::result<T>
tag_invoke( try_value_to_tag<T>, const value& jv, const use_allocator_t<Alloc>& ctx, const FullContext& full_ctx )
{

    array const* arr = jv.if_array();
    if( !arr )
        return {
            boost::system::in_place_error,
            make_error_code(boost::system::errc::invalid_argument) };

    T result(ctx.allocator);
    auto ins = std::inserter(result, result.end());
    for( value const& val: *arr )
    {
        using ValueType = typename T::value_type;
        auto elem_res = try_value_to<ValueType>( val, full_ctx );
        if( elem_res.has_error() )
            return {boost::system::in_place_error, elem_res.error()};
        *ins++ = std::move(*elem_res);
    }
    return result;
}

int
main(int, char**)
{

    value const jv = { 1, 2, 3, 4, 5, 6, 7, 8 };

    unsigned char buf[1024];
    pmr::monotonic_buffer_resource mr( buf, sizeof(buf) );

    auto v = value_to< pmr::vector<int> >(
        jv,
        use_allocator( pmr::polymorphic_allocator<int>(&mr) ) );

    assert( v.size() == jv.as_array().size() );

    for( auto i = 0u; i < v.size(); ++i )
        assert( v[i] == jv.at(i).to_number<int>() );

    return EXIT_SUCCESS;
}

//]
