// Copyright 2018-2024 Emil Dotchevski and Reverge Studios, Inc.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef BOOST_LEAF_TEST_SINGLE_HEADER
#   include "leaf.hpp"
#else
#   include <boost/leaf/result.hpp>
#   include <boost/leaf/handle_errors.hpp>
#endif

#include "lightweight_test.hpp"

namespace leaf = boost::leaf;

static_assert(std::is_same<void, leaf::result<void>::value_type>::value, "Bad value_type");
static_assert(std::is_same<int, leaf::result<int>::value_type>::value, "Bad value_type");
static_assert(std::is_same<int const, leaf::result<int const>::value_type>::value, "Bad value_type");
static_assert(std::is_same<int &, leaf::result<int &>::value_type>::value, "Bad value_type");
static_assert(std::is_same<int const &, leaf::result<int const &>::value_type>::value, "Bad value_type");

struct val
{
    static int id_count;
    static int count;
    int id;
    float a = 0;
    int b = 0;

    val( float a, int b ) noexcept:
        id(++id_count),
        a(a),
        b(b)
    {
        ++count;
    }

    val() noexcept:
        id(++id_count)
    {
        ++count;
    }

    val( val const & x ) noexcept:
        id(x.id)
    {
        ++count;
    }

    val( val && x ) noexcept:
        id(x.id)
    {
        ++count;
    }

    ~val() noexcept
    {
        --count;
    }

    friend bool operator==( val const & a, val const & b ) noexcept
    {
        return a.id == b.id;
    }

    friend std::ostream & operator<<( std::ostream & os, val const & v ) noexcept
    {
        return os << v.id;
    }
};
int val::count = 0;
int val::id_count = 0;

struct err
{
    static int count;

    err()
    {
        ++count;
    }

    err( err const & )
    {
        ++count;
    }

    err( err && )
    {
        ++count;
    }

    ~err()
    {
        --count;
    }
};
int err::count = 0;
struct e_err { err value; };

bool eq_value( leaf::result<val> & r, val v )
{
    leaf::result<val> const & cr = r;
    val const & cv = v;
    return
        r.value() == v &&
        cr.value() == cv &&
        *r.operator->() == v &&
        *cr.operator->() == cv &&
        *r == v &&
        *cr == cv;
}

int main()
{
    { // value default -> move
        leaf::result<val> r1;
        BOOST_TEST(r1);
        BOOST_TEST_EQ(err::count, 0);
        BOOST_TEST_EQ(val::count, 1);
        leaf::result<val> r2 = std::move(r1);
        BOOST_TEST(r2);
        BOOST_TEST_EQ(err::count, 0);
        BOOST_TEST_EQ(val::count, 2);
    }
    BOOST_TEST_EQ(err::count, 0);
    BOOST_TEST_EQ(val::count, 0);
    { // value move -> move
        leaf::result<val> r1 = val();
        BOOST_TEST(r1);
        BOOST_TEST_EQ(err::count, 0);
        BOOST_TEST_EQ(val::count, 1);
        leaf::result<val> r2 = std::move(r1);
        BOOST_TEST(r2);
        BOOST_TEST_EQ(err::count, 0);
        BOOST_TEST_EQ(val::count, 2);
    }
    BOOST_TEST_EQ(err::count, 0);
    BOOST_TEST_EQ(val::count, 0);
    { // value copy -> move
        val v;
        leaf::result<val> r1 = v;
        BOOST_TEST(r1);
        BOOST_TEST_EQ(err::count, 0);
        BOOST_TEST_EQ(val::count, 2);
        BOOST_TEST(eq_value(r1, v));
        leaf::result<val> r2 = std::move(r1);
        BOOST_TEST(r2);
        BOOST_TEST_EQ(err::count, 0);
        BOOST_TEST_EQ(val::count, 3);
        BOOST_TEST(eq_value(r2, v));
    }
    BOOST_TEST_EQ(err::count, 0);
    BOOST_TEST_EQ(val::count, 0);
    { // value emplace -> move
        leaf::result<val> r1 = { 42.0f, 42 };
        BOOST_TEST(r1);
        BOOST_TEST_EQ(r1.value().a, 42.0f);
        BOOST_TEST_EQ(r1.value().b, 42);
        BOOST_TEST_EQ(err::count, 0);
        BOOST_TEST_EQ(val::count, 1);
        leaf::result<val> r2 = std::move(r1);
        BOOST_TEST(r2);
        BOOST_TEST_EQ(err::count, 0);
        BOOST_TEST_EQ(val::count, 2);
    }
    BOOST_TEST_EQ(err::count, 0);
    BOOST_TEST_EQ(val::count, 0);

    { // value default -> assign-move
        leaf::result<val> r1;
        BOOST_TEST(r1);
        BOOST_TEST_EQ(err::count, 0);
        BOOST_TEST_EQ(val::count, 1);
        leaf::result<val> r2; r2=std::move(r1);
        BOOST_TEST(r2);
        BOOST_TEST_EQ(err::count, 0);
        BOOST_TEST_EQ(val::count, 2);
    }
    BOOST_TEST_EQ(err::count, 0);
    BOOST_TEST_EQ(val::count, 0);
    { // value move -> assign-move
        leaf::result<val> r1 = val();
        BOOST_TEST(r1);
        BOOST_TEST_EQ(err::count, 0);
        BOOST_TEST_EQ(val::count, 1);
        leaf::result<val> r2; r2=std::move(r1);
        BOOST_TEST(r2);
        BOOST_TEST_EQ(err::count, 0);
        BOOST_TEST_EQ(val::count, 2);
    }
    BOOST_TEST_EQ(err::count, 0);
    BOOST_TEST_EQ(val::count, 0);
    { // value copy -> assign-move
        val v;
        leaf::result<val> r1 = v;
        BOOST_TEST(r1);
        BOOST_TEST_EQ(err::count, 0);
        BOOST_TEST_EQ(val::count, 2);
        BOOST_TEST(eq_value(r1, v));
        leaf::result<val> r2; r2=std::move(r1);
        BOOST_TEST(r2);
        BOOST_TEST_EQ(err::count, 0);
        BOOST_TEST_EQ(val::count, 3);
        BOOST_TEST(eq_value(r2, v));
    }
    BOOST_TEST_EQ(err::count, 0);
    BOOST_TEST_EQ(val::count, 0);
    { // value emplace -> assign-move
        leaf::result<val> r1 = { 42.0f, 42 };
        BOOST_TEST(r1);
        BOOST_TEST_EQ(err::count, 0);
        BOOST_TEST_EQ(val::count, 1);
        leaf::result<val> r2; r2=std::move(r1);
        BOOST_TEST(r2);
        BOOST_TEST_EQ(err::count, 0);
        BOOST_TEST_EQ(val::count, 2);
    }
    BOOST_TEST_EQ(err::count, 0);
    BOOST_TEST_EQ(val::count, 0);

    // ^^ value ^^
    // vv error vv

    { // error move -> move
        leaf::context<e_err>  ctx;
        auto active_context = activate_context(ctx);
        leaf::result<val> r1 = leaf::new_error( e_err { } );
        BOOST_TEST(!r1);
        BOOST_TEST_EQ(err::count, 1);
        BOOST_TEST_EQ(val::count, 0);
        leaf::error_id r1e = r1.error();
        leaf::result<val> r2 = std::move(r1);
        leaf::error_id r2e = r2.error();
        BOOST_TEST_EQ(r1e, r2e);
        BOOST_TEST(!r2);
    }
    BOOST_TEST_EQ(err::count, 0);
    BOOST_TEST_EQ(val::count, 0);
    { // error copy -> move
        leaf::context<e_err>  ctx;
        auto active_context = activate_context(ctx);
        leaf::error_id err = leaf::new_error( e_err{ } );
        leaf::result<val> r1 = err;
        BOOST_TEST(!r1);
        BOOST_TEST_EQ(err::count, 1);
        BOOST_TEST_EQ(val::count, 0);
        leaf::error_id r1e = r1.error();
        leaf::result<val> r2 = std::move(r1);
        leaf::error_id r2e = r2.error();
        BOOST_TEST_EQ(r1e, r2e);
        BOOST_TEST(!r2);
    }
    BOOST_TEST_EQ(err::count, 0);
    BOOST_TEST_EQ(val::count, 0);

    { // error move -> assign move
        leaf::context<e_err>  ctx;
        ctx.activate();
        leaf::result<val> r1 = leaf::new_error( e_err { } );
        ctx.deactivate();
        BOOST_TEST(!r1);
        BOOST_TEST_EQ(err::count, 1);
        BOOST_TEST_EQ(val::count, 0);
        leaf::error_id r1e = r1.error();
        leaf::result<val> r2; r2=std::move(r1);
        leaf::error_id r2e = r2.error();
        BOOST_TEST_EQ(r1e, r2e);
        BOOST_TEST(!r2);
        {
            val x;
            BOOST_TEST(ctx.handle_error<val>(r2.error(), [&]{ return x; }) == x);
        }
        BOOST_TEST_EQ(err::count, 1);
        BOOST_TEST_EQ(val::count, 0);
    }
    BOOST_TEST_EQ(err::count, 0);
    BOOST_TEST_EQ(val::count, 0);
    { // error copy -> assign move
        leaf::context<e_err>  ctx;
        auto active_context = activate_context(ctx);
        leaf::error_id err = leaf::new_error( e_err{ } );
        leaf::result<val> r1 = err;
        BOOST_TEST(!r1);
        BOOST_TEST_EQ(err::count, 1);
        BOOST_TEST_EQ(val::count, 0);
        leaf::error_id r1e = r1.error();
        leaf::result<val> r2; r2=std::move(r1);
        leaf::error_id r2e = r2.error();
        BOOST_TEST_EQ(r1e, r2e);
        BOOST_TEST(!r2);
    }
    BOOST_TEST_EQ(err::count, 0);
    BOOST_TEST_EQ(val::count, 0);

#if BOOST_LEAF_CFG_CAPTURE
    { // error move -> capture -> move
        leaf::result<val> r1 = leaf::try_capture_all(
            []()->leaf::result<val>
            {
                return leaf::new_error(e_err{ });
            });
        BOOST_TEST(!r1);
        BOOST_TEST_EQ(err::count, 1);
        BOOST_TEST_EQ(val::count, 0);
        leaf::error_id r1e = r1.error();
        leaf::result<val> r2 = std::move(r1);
        leaf::error_id r2e = r2.error();
        BOOST_TEST_EQ(r1e, r2e);
        BOOST_TEST(!r2);
        BOOST_TEST(!r1);
    }
    BOOST_TEST_EQ(err::count, 0);
    BOOST_TEST_EQ(val::count, 0);
    { // error copy -> capture -> move
        leaf::result<val> r1 = leaf::try_capture_all(
            []()->leaf::result<val>
            {
                leaf::error_id err = leaf::new_error( e_err{ } );
                return leaf::result<val>(err);
            });
        BOOST_TEST(!r1);
        BOOST_TEST_EQ(err::count, 1);
        BOOST_TEST_EQ(val::count, 0);
        leaf::error_id r1e = r1.error();
        leaf::result<val> r2 = std::move(r1);
        leaf::error_id r2e = r2.error();
        BOOST_TEST_EQ(r1e, r2e);
        BOOST_TEST(!r2);
        BOOST_TEST(!r1);
    }
    BOOST_TEST_EQ(err::count, 0);
    BOOST_TEST_EQ(val::count, 0);

    { // error move -> capture -> assign-move
        leaf::result<val> r1 = leaf::try_capture_all(
            []()->leaf::result<val>
            {
                return leaf::new_error(e_err{ });
            });
        BOOST_TEST(!r1);
        BOOST_TEST_EQ(err::count, 1);
        BOOST_TEST_EQ(val::count, 0);
        leaf::error_id r1e = r1.error();
        leaf::result<val> r2; r2=std::move(r1);
        leaf::error_id r2e = r2.error();
        BOOST_TEST_EQ(r1e, r2e);
        BOOST_TEST(!r2);
        BOOST_TEST(!r1);
    }
    BOOST_TEST_EQ(err::count, 0);
    BOOST_TEST_EQ(val::count, 0);
    { // error copy -> capture -> assign-move
        leaf::result<val> r1 = leaf::try_capture_all(
            []()->leaf::result<val>
            {
                leaf::error_id err = leaf::new_error( e_err{ } );
                return leaf::result<val>(err);
            });
        BOOST_TEST(!r1);
        BOOST_TEST_EQ(err::count, 1);
        BOOST_TEST_EQ(val::count, 0);
        leaf::error_id r1e = r1.error();
        leaf::result<val> r2; r2=std::move(r1);
        leaf::error_id r2e = r2.error();
        BOOST_TEST_EQ(r1e, r2e);
        BOOST_TEST(!r2);
        BOOST_TEST(!r1);
    }
    BOOST_TEST_EQ(err::count, 0);
    BOOST_TEST_EQ(val::count, 0);
#endif

    // ^^ result<T> ^^

    ////////////////////////////////////////

    // vv result<void> vv

    { // void default -> move
        leaf::result<void> r1;
        BOOST_TEST(r1);
        r1.value();
        BOOST_TEST(r1.operator->());
        *r1;
        leaf::result<void> r2 = std::move(r1);
        BOOST_TEST(r2);
        r2.value();
        BOOST_TEST(r2.operator->());
        *r2;
    }

    { // void default -> assign-move
        leaf::result<void> r1;
        BOOST_TEST(r1);
        r1.value();
        BOOST_TEST(r1.operator->());
        *r1;
        leaf::result<void> r2; r2=std::move(r1);
        BOOST_TEST(r2);
        r2.value();
        BOOST_TEST(r2.operator->());
        *r2;
    }

    // ^^ void default ^^
    // vv void error vv

    { // void error move -> move
        leaf::context<e_err>  ctx;
        auto active_context = activate_context(ctx);
        leaf::result<void> r1 = leaf::new_error( e_err { } );
        BOOST_TEST(!r1);
        BOOST_TEST(!r1.operator->());
        BOOST_TEST_EQ(err::count, 1);
        leaf::error_id r1e = r1.error();
        leaf::result<void> r2 = std::move(r1);
        leaf::error_id r2e = r2.error();
        BOOST_TEST_EQ(r1e, r2e);
        BOOST_TEST(!r2);
        BOOST_TEST(r2.operator->() == 0);
    }
    BOOST_TEST_EQ(err::count, 0);
    { // void error copy -> move
        leaf::context<e_err>  ctx;
        auto active_context = activate_context(ctx);
        leaf::error_id err = leaf::new_error( e_err{ } );
        leaf::result<void> r1 = err;
        BOOST_TEST(!r1);
        BOOST_TEST(!r1.operator->());
        BOOST_TEST_EQ(err::count, 1);
        leaf::error_id r1e = r1.error();
        leaf::result<void> r2 = std::move(r1);
        leaf::error_id r2e = r2.error();
        BOOST_TEST_EQ(r1e, r2e);
        BOOST_TEST(!r2);
        BOOST_TEST(r2.operator->() == 0);
    }
    BOOST_TEST_EQ(err::count, 0);

    { // void error move -> assign move
        leaf::context<e_err>  ctx;
        ctx.activate();
        leaf::result<void> r1 = leaf::new_error( e_err { } );
        ctx.deactivate();
        BOOST_TEST(!r1);
        BOOST_TEST(!r1.operator->());
        BOOST_TEST_EQ(err::count, 1);
        leaf::error_id r1e = r1.error();
        leaf::result<void> r2; r2=std::move(r1);
        leaf::error_id r2e = r2.error();
        BOOST_TEST_EQ(r1e, r2e);
        BOOST_TEST(!r2);
        BOOST_TEST(r2.operator->() == 0);
        ctx.handle_error<void>(r2.error(), []{ });
        BOOST_TEST_EQ(err::count, 1);
    }
    BOOST_TEST_EQ(err::count, 0);
    { // void error copy -> assign move
        leaf::context<e_err>  ctx;
        auto active_context = activate_context(ctx);
        leaf::error_id err = leaf::new_error( e_err{ } );
        leaf::result<void> r1 = err;
        BOOST_TEST(!r1);
        BOOST_TEST(!r1.operator->());
        BOOST_TEST_EQ(err::count, 1);
        leaf::error_id r1e = r1.error();
        leaf::result<void> r2; r2=std::move(r1);
        leaf::error_id r2e = r2.error();
        BOOST_TEST_EQ(r1e, r2e);
        BOOST_TEST(!r2);
        BOOST_TEST(r2.operator->() == 0);
    }
    BOOST_TEST_EQ(err::count, 0);

#if BOOST_LEAF_CFG_CAPTURE
    { // void error move -> capture -> move
        leaf::result<void> r1 = leaf::try_capture_all(
            []()->leaf::result<void>
            {
                return leaf::new_error(e_err{ });
            });
        BOOST_TEST(!r1);
        BOOST_TEST(!r1.operator->());
        BOOST_TEST_EQ(err::count, 1);
        leaf::error_id r1e = r1.error();
        leaf::result<void> r2 = std::move(r1);
        leaf::error_id r2e = r2.error();
        BOOST_TEST_EQ(r1e, r2e);
        BOOST_TEST(!r2);
        BOOST_TEST(r2.operator->() == 0);
    }
    BOOST_TEST_EQ(err::count, 0);
    { // void error copy -> capture -> move
        leaf::result<void> r1 = leaf::try_capture_all(
            []()->leaf::result<void>
            {
                leaf::error_id err = leaf::new_error( e_err{ } );
                return leaf::result<void>(err);
            });
        BOOST_TEST(!r1);
        BOOST_TEST(!r1.operator->());
        BOOST_TEST_EQ(err::count, 1);
        leaf::error_id r1e = r1.error();
        leaf::result<void> r2 = std::move(r1);
        leaf::error_id r2e = r2.error();
        BOOST_TEST_EQ(r1e, r2e);
        BOOST_TEST(!r2);
        BOOST_TEST(r2.operator->() == 0);
    }
    BOOST_TEST_EQ(err::count, 0);

    { // void error move -> capture -> assign-move
        leaf::result<void> r1 = leaf::try_capture_all(
            []()->leaf::result<void>
            {
                return leaf::new_error(e_err{ });
            });
        BOOST_TEST(!r1);
        BOOST_TEST(!r1.operator->());
        BOOST_TEST_EQ(err::count, 1);
        leaf::error_id r1e = r1.error();
        leaf::result<void> r2; r2=std::move(r1);
        leaf::error_id r2e = r2.error();
        BOOST_TEST_EQ(r1e, r2e);
        BOOST_TEST(!r2);
        BOOST_TEST(r2.operator->() == 0);
    }
    BOOST_TEST_EQ(err::count, 0);
    { // void error copy -> capture -> assign-move
        leaf::result<void> r1 = leaf::try_capture_all(
            []()->leaf::result<void>
            {
                leaf::error_id err = leaf::new_error( e_err{ } );
                return leaf::result<void>(err);
            });
        BOOST_TEST(!r1);
        BOOST_TEST(!r1.operator->());
        BOOST_TEST_EQ(err::count, 1);
        leaf::error_id r1e = r1.error();
        leaf::result<void> r2; r2=std::move(r1);
        leaf::error_id r2e = r2.error();
        BOOST_TEST_EQ(r1e, r2e);
        BOOST_TEST(!r2);
        BOOST_TEST(r2.operator->() == 0);
    }
    BOOST_TEST_EQ(err::count, 0);
#endif

    {
        leaf::result<int> r = leaf::error_id();
        BOOST_TEST(!r);
        BOOST_TEST_EQ(val::count, 0);
        BOOST_TEST_EQ(err::count, 0);
    }
    BOOST_TEST_EQ(val::count, 0);
    BOOST_TEST_EQ(err::count, 0);

    {
        leaf::result<void> r = leaf::error_id();
        BOOST_TEST(!r);
        BOOST_TEST_EQ(val::count, 0);
        BOOST_TEST_EQ(err::count, 0);
    }
    BOOST_TEST_EQ(val::count, 0);
    BOOST_TEST_EQ(err::count, 0);

    {
        leaf::result<void> r;
        BOOST_TEST(r);
        leaf::result<val> r1 = r.error();
        BOOST_TEST_EQ(val::count, 0);
        BOOST_TEST(!r1);
        leaf::error_id id = r.error();
        BOOST_TEST(!id);
    }
    BOOST_TEST_EQ(val::count, 0);

    {
        leaf::result<val> r;
        BOOST_TEST(r);
        leaf::result<void> r1 = r.error();
        BOOST_TEST(!r1);
        leaf::error_id id = r.error();
        BOOST_TEST(!id);
        BOOST_TEST_EQ(val::count, 1);
    }
    BOOST_TEST_EQ(val::count, 0);

    {
        leaf::result<val> r;
        BOOST_TEST(r);
        leaf::result<float> r1 = r.error();
        BOOST_TEST(!r1);
        leaf::error_id id = r.error();
        BOOST_TEST(!id);
        BOOST_TEST_EQ(val::count, 1);
    }
    BOOST_TEST_EQ(val::count, 0);

#if BOOST_LEAF_CFG_STD_STRING
    { // Initialization forwarding constructor
        leaf::result<std::string> r = "hello";
        BOOST_TEST(r);
        BOOST_TEST_EQ(r.value(), "hello");
    }
#endif

#if BOOST_LEAF_CFG_STD_STRING
    { // Initialization forwarding constructor
        leaf::result<std::string> r; r = "hello";
        BOOST_TEST(r);
        BOOST_TEST_EQ(r.value(), "hello");
    }
#endif

    return boost::report_errors();
}
