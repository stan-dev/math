//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

#ifndef BOOST_URL_EXTRA_TEST_SUITE_HPP
#define BOOST_URL_EXTRA_TEST_SUITE_HPP

#if defined(_MSC_VER)
# pragma once
#endif

#include <boost/current_function.hpp>
#include <cctype>
#include <sstream>
#include <type_traits>

//  This is a derivative work
//  Copyright 2002-2018 Peter Dimov
//  Copyright (c) 2002, 2009, 2014 Peter Dimov
//  Copyright (2) Beman Dawes 2010, 2011
//  Copyright (3) Ion Gaztanaga 2013
//
//  Copyright 2018 Glen Joseph Fernandes
//  (glenjofe@gmail.com)

namespace test_suite {

//------------------------------------------------

struct any_suite
{
    virtual ~any_suite() = 0;
    virtual char const* name() const noexcept = 0;
    virtual void run() const = 0;
};

//------------------------------------------------

struct suites
{
    virtual ~suites() = default;

    using iterator = any_suite const* const*;
    virtual void insert(any_suite const&) = 0;
    virtual iterator begin() const noexcept = 0;
    virtual iterator end() const noexcept = 0;

    // DEPRECATED
    virtual void sort() = 0;

    static suites& instance() noexcept;
};

//------------------------------------------------

template<class T>
class suite : public any_suite
{
    char const* name_;

public:
    explicit
    suite(char const* name) noexcept
        : name_(name)
    {
        suites::instance().insert(*this);
    }

    char const*
    name() const noexcept override
    {
        return name_;
    }

    void
    run() const override
    {
        T().run();
    }
};

//------------------------------------------------

class any_runner
{
    any_runner* prev_;

    static any_runner*&
    instance_impl() noexcept;

public:
    static any_runner& instance() noexcept;

    any_runner() noexcept;
    virtual ~any_runner();

    virtual void run(any_suite const& test) = 0;
    virtual void note(char const* msg) = 0;
    virtual bool test(bool cond,
        char const* expr, char const* func,
        char const* file, int line) = 0;
    virtual std::ostream& log() noexcept = 0;
};

//------------------------------------------------

namespace detail {

// In the comparisons below, it is possible that
// T and U are signed and unsigned integer types,
// which generates warnings in some compilers.
// A cleaner fix would require common_type trait
// or some meta-programming, which would introduce
// a dependency on Boost.TypeTraits. To avoid
// the dependency we just disable the warnings.
#if defined(__clang__) && defined(__has_warning)
# if __has_warning("-Wsign-compare")
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wsign-compare"
# endif
#elif defined(_MSC_VER)
# pragma warning(push)
# pragma warning(disable: 4389)
# pragma warning(disable: 4018)
#elif defined(__GNUC__) && !(defined(__INTEL_COMPILER) || defined(__ICL) || defined(__ICC) || defined(__ECC)) && (__GNUC__ * 100 + __GNUC_MINOR__) >= 406
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wsign-compare"
#endif

template<class...>
struct make_void
{
    using type = void;
};

template<class... Ts>
using void_t = typename make_void<Ts...>::type;

template <class T, class = void>
struct is_streamable : std::false_type
{};

template <class T>
struct is_streamable<
    T, void_t<decltype(std::declval<
        std::ostream&>() << std::declval<T&>())
    > > : std::true_type
{};

template <class T>
auto
test_output_impl(T const& v) ->
    typename std::enable_if<
        is_streamable<T>::value,
        T const&>::type
{
    return v;
}

template <class T>
auto
test_output_impl(T const&) ->
    typename std::enable_if<
        ! is_streamable<T>::value,
        std::string>::type
{
    return "?";
}

// specialize test output for char pointers to avoid printing as cstring
template<class T>
       const void* test_output_impl(T volatile* v) { return const_cast<T*>(v); }
inline const void* test_output_impl(const char* v) { return v; }
inline const void* test_output_impl(const unsigned char* v) { return v; }
inline const void* test_output_impl(const signed char* v) { return v; }
inline const void* test_output_impl(char* v) { return v; }
inline const void* test_output_impl(unsigned char* v) { return v; }
inline const void* test_output_impl(signed char* v) { return v; }
inline const void* test_output_impl(std::nullptr_t) { return nullptr; }

// print chars as numeric
inline int test_output_impl( signed char const& v ) { return v; }
inline unsigned test_output_impl( unsigned char const& v ) { return v; }

// Whether wchar_t is signed is implementation-defined
template<bool Signed> struct lwt_long_type {};
template<> struct lwt_long_type<true> { typedef long type; };
template<> struct lwt_long_type<false> { typedef unsigned long type; };
inline lwt_long_type<
    (static_cast<wchar_t>(-1) < static_cast<wchar_t>(0))
        >::type test_output_impl( wchar_t const& v ) { return v; }

#if !defined( BOOST_NO_CXX11_CHAR16_T )
inline unsigned long test_output_impl( char16_t const& v ) { return v; }
#endif

#if !defined( BOOST_NO_CXX11_CHAR32_T )
inline unsigned long test_output_impl( char32_t const& v ) { return v; }
#endif

#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable: 4996)
#endif

inline std::string test_output_impl( char const& v )
{
    if( std::isprint( static_cast<unsigned char>( v ) ) )
    {
        return { 1, v };
    }
    else
    {
        char buffer[ 8 ];
        std::snprintf( buffer, sizeof(buffer), "\\x%02X", static_cast<unsigned char>( v ) );
        return buffer;
    }
}

bool
test_impl(
    bool cond,
    char const* expr,
    char const* func,
    char const* file,
    int line);

void
throw_failed_impl(
    const char* expr,
    char const* excep,
    char const* func,
    char const* file,
    int line);

void
no_throw_failed_impl(
    const char* expr,
    char const* excep,
    char const* func,
    char const* file,
    int line);

void
no_throw_failed_impl(
    const char* expr,
    char const* func,
    char const* file,
    int line);

struct lw_test_eq
{
    template <typename T, typename U>
    bool operator()(const T& t, const U& u) const
    {
        return t == u;
    }
};

struct lw_test_ne
{

    template <typename T, typename U>
    bool operator()(const T& t, const U& u) const
    {
        return t != u;
    }
};

struct lw_test_lt
{
    template <typename T, typename U>
    bool operator()(const T& t, const U& u) const
    {
        return t < u;
    }
};

struct lw_test_gt
{
    template <typename T, typename U>
    bool operator()(const T& t, const U& u) const
    {
        return t > u;
    }
};

struct lw_test_le
{
    template <typename T, typename U>
    bool operator()(const T& t, const U& u) const
    {
        return t <= u;
    }
};

struct lw_test_ge
{
    template <typename T, typename U>
    bool operator()(const T& t, const U& u) const
    {
        return t >= u;
    }
};

// lwt_predicate_name

template<class T> char const * lwt_predicate_name( T const& )
{
    return "~=";
}

inline char const * lwt_predicate_name( lw_test_eq const& )
{
    return "==";
}

inline char const * lwt_predicate_name( lw_test_ne const& )
{
    return "!=";
}

inline char const * lwt_predicate_name( lw_test_lt const& )
{
    return "<";
}

inline char const * lwt_predicate_name( lw_test_le const& )
{
    return "<=";
}

inline char const * lwt_predicate_name( lw_test_gt const& )
{
    return ">";
}

inline char const * lwt_predicate_name( lw_test_ge const& )
{
    return ">=";
}

//------------------------------------------------

template<class Pred, class T, class U>
bool
test_with_impl(
    Pred pred,
    char const* expr1,
    char const* expr2,
    char const* func,
    char const* file,
    int line,
    T const& t, U const& u)
{
    if(pred(t, u))
    {
        any_runner::instance().test(
            true, "", func, file, line);
        return true;
    }
    std::stringstream ss;
    ss <<
        "\"" << test_output_impl(t) << "\" " <<
        lwt_predicate_name(pred) <<
        " \"" << test_output_impl(u) << "\" (" <<
        expr1 << " " <<
        lwt_predicate_name(pred) <<
        " " << expr2 << ")";
    any_runner::instance().test(
        false, ss.str().c_str(), func, file, line);
    return false;
}

#if defined(__clang__) && defined(__has_warning)
# if __has_warning("-Wsign-compare")
#  pragma clang diagnostic pop
# endif
#elif defined(_MSC_VER)
# pragma warning(pop)
#elif defined(__GNUC__) && !(defined(__INTEL_COMPILER) || defined(__ICL) || defined(__ICC) || defined(__ECC)) && (__GNUC__ * 100 + __GNUC_MINOR__) >= 406
# pragma GCC diagnostic pop
#endif

//------------------------------------------------

struct log_type
{
    template<class T>
    friend
    std::ostream&
    operator<<(
        log_type const&, T&& t)
    {
        return any_runner::instance().log() << t;
    }
};

//------------------------------------------------

} // detail

/** Log output to the current suite
*/
constexpr detail::log_type log{};

#define BOOST_TEST(expr) ( \
    ::test_suite::detail::test_impl( \
        (expr) ? true : false, #expr, \
            BOOST_CURRENT_FUNCTION, __FILE__, __LINE__ ) )

#define BOOST_ERROR(msg) \
    ::test_suite::detail::test_impl( \
        false, msg, BOOST_CURRENT_FUNCTION, __FILE__, __LINE__ )

#define BOOST_TEST_WITH(expr1,expr2,predicate) ( \
    ::test_suite::detail::test_with_impl( \
        predicate, #expr1, #expr2, BOOST_CURRENT_FUNCTION, \
        __FILE__, __LINE__, expr1, expr2) )

#define BOOST_TEST_EQ(expr1,expr2) \
    BOOST_TEST_WITH( expr1, expr2, \
        ::test_suite::detail::lw_test_eq() )

#define BOOST_TEST_CSTR_EQ(expr1,expr2) \
    BOOST_TEST_EQ( string_view(expr1), string_view(expr2) )

#define BOOST_TEST_NE(expr1,expr2) \
    BOOST_TEST_WITH( expr1, expr2, \
        ::test_suite::detail::lw_test_ne() )

#define BOOST_TEST_LT(expr1,expr2) \
    BOOST_TEST_WITH( expr1, expr2, \
        ::test_suite::detail::lw_test_lt() )

#define BOOST_TEST_LE(expr1,expr2) \
    BOOST_TEST_WITH( expr1, expr2, \
        ::test_suite::detail::lw_test_le() )

#define BOOST_TEST_GT(expr1,expr2) \
    BOOST_TEST_WITH( expr1, expr2, \
        ::test_suite::detail::lw_test_gt() )

#define BOOST_TEST_GE(expr1,expr2) \
    BOOST_TEST_WITH( expr1, expr2, \
        ::test_suite::detail::lw_test_ge() )

#define BOOST_TEST_PASS() BOOST_TEST(true)

#define BOOST_TEST_FAIL() BOOST_TEST(false)

#define BOOST_TEST_NOT(expr) BOOST_TEST(!(expr))

#ifndef BOOST_NO_EXCEPTIONS
# define BOOST_TEST_THROWS( expr, ex ) \
    try { \
        (void)(expr); \
        ::test_suite::detail::throw_failed_impl( \
            #expr, #ex, BOOST_CURRENT_FUNCTION, \
            __FILE__, __LINE__); \
    } \
    catch(ex const&) { \
        BOOST_TEST_PASS(); \
    } \
    catch(...) { \
        ::test_suite::detail::throw_failed_impl( \
            #expr, #ex, BOOST_CURRENT_FUNCTION, \
            __FILE__, __LINE__); \
    }
   //
#else
   #define BOOST_TEST_THROWS( expr, ex )
#endif

#ifndef BOOST_NO_EXCEPTIONS
# define BOOST_TEST_NO_THROW( expr ) \
    try { \
        (void)(expr); \
        BOOST_TEST_PASS(); \
    } catch (const std::exception& e) { \
        ::test_suite::detail::no_throw_failed_impl( \
            #expr, e.what(), BOOST_CURRENT_FUNCTION, \
            __FILE__, __LINE__); \
    } catch (...) { \
        ::test_suite::detail::no_throw_failed_impl( \
            #expr, BOOST_CURRENT_FUNCTION, \
            __FILE__, __LINE__); \
    }
    //
#else
# define BOOST_TEST_NO_THROW( expr ) ( [&]{ expr; return true; }() )
#endif

#define TEST_SUITE(type, name) \
    static ::test_suite::suite<type> type##_(name)

} // test_suite

#endif
