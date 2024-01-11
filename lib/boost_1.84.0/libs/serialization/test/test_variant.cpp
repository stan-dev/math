/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// test_variant.cpp
// test of non-intrusive serialization of variant types
//
// copyright (c) 2005
// troy d. straszheim <troy@resophonic.com>
// http://www.resophonic.com
//
// copyright (c) 2023
// Robert Ramey <ramey@rrsd.com>
// http://www.rrsd.com
//
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org for updates, documentation, and revision history.
//
// thanks to Robert Ramey and Peter Dimov.
//

#include <cstddef> // NULL
#include <cstdio> // remove
#include <fstream>

#include <boost/config.hpp>
#if BOOST_CXX_VERSION > 199711L // only include floating point if C++ version >= C++11
#include <boost/math/special_functions/next.hpp>
#endif

#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std {
    using ::remove;
}
#endif

#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/identity.hpp>

#if defined(_MSC_VER) && (_MSC_VER <= 1020)
#  pragma warning (disable : 4786) // too long name, harmless warning
#endif

#include "test_tools.hpp"

#include <boost/archive/archive_exception.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/throw_exception.hpp>

#include <boost/variant/static_visitor.hpp>

namespace boost {
    template<typename ResultType> class static_visitor;
}

#include "A.hpp"
#include "A.ipp"

class are_equal
    : public boost::static_visitor<bool>
{
public:
    typedef bool result_type;
    // note extra rigamarole for compilers which don't support
    // partial function template ordering - specifically msvc 6.x
    struct same {
        template<class T, class U>
        static bool invoke(const T & t, const U & u){
            return t == u;
        }
    };

    struct not_same {
        template<class T, class U>
        static bool invoke(const T &, const U &){
            return false;
        }
    };

    template <class T, class U>
    bool operator()( const T & t, const U & u) const
    {
        typedef typename boost::mpl::eval_if<boost::is_same<T, U>,
            boost::mpl::identity<same>,
            boost::mpl::identity<not_same>
        >::type type;
        return type::invoke(t, u);
    }

    template <class T, class U>
    bool operator()(T * const & t,  U * const & u) const
    {
        return this->operator()(*t, *u);
    }

    bool operator()( const float & lhs, const float & rhs ) const
    {
        #if BOOST_CXX_VERSION > 199711L // only include floating point if C++ version >= C++11
        return std::abs( boost::math::float_distance(lhs, rhs) ) < 2;
        #else
        return true;
        #endif
    }
    bool operator()( const double & lhs, const double & rhs ) const
    {
        #if BOOST_CXX_VERSION > 199711L // only include floating point if C++ version >= C++11
        return std::abs( boost::math::float_distance(lhs, rhs) ) < 2;
        #else
        return true;
        #endif
    }
};

template<class Variant>
bool test_type(const Variant & v){
    const char * testfile = boost::archive::tmpnam(NULL);
    BOOST_REQUIRE(testfile != NULL);
    {
        test_ostream os(testfile, TEST_STREAM_FLAGS);
        test_oarchive oa(os, TEST_ARCHIVE_FLAGS);
        oa << boost::serialization::make_nvp("written", v);
    }

    Variant vx;
    {
        test_istream is(testfile, TEST_STREAM_FLAGS);
        test_iarchive ia(is, TEST_ARCHIVE_FLAGS);
        BOOST_TRY {
            ia >> boost::serialization::make_nvp("written", vx);
            BOOST_CHECK(visit(are_equal(), v, vx));
        }
        BOOST_CATCH(boost::archive::archive_exception const& e) {
            return false;
        }
        BOOST_CATCH_END
    }
    std::remove(testfile);
    return true;
}

template<class Variant>
void test(Variant & v)
{
    // uninitialized
    test_type(v);
    v = false;
    test_type(v);
    v = 1;
    test_type(v);
    v = (float) 2.3;
    test_type(v);
    v = (double) 6.4;
    test_type(v);
    v = A();
    test_type(v);
    v = std::string("we can't stop here, this is Bat Country");
    test_type(v);
}

#include <boost/serialization/variant.hpp>

int test_main( int /* argc */, char* /* argv */[] ){

    // boost::variant - compatible with C++03
    {
        boost::variant<bool, int, float, double, A, std::string> v;
        test(v);
        const A a;
        boost::variant<bool, int, float, double, const A *, std::string> v1 = & a;
        test_type(v1);
    }

    // boost::variant2/variant requires C++ 11
    #if BOOST_CXX_VERSION >= 201103L
    {
        boost::variant2::variant<bool, int, float, double, A, std::string> v;
        test(v);
        const A a;
        boost::variant2::variant<bool, int, float, double, const A *, std::string> v1 = & a;
        test_type(v1);
    }
    #endif

    // std::variant reqires C++ 17 or more
    #ifndef BOOST_NO_CXX17_HDR_VARIANT
    {
        std::variant<bool, int, float, double, A, std::string> v;
        test(v);
        const A a;
        std::variant<bool, int, float, double, const A *, std::string> v1 = & a;
        test_type(v1);
    }
    #endif

    return EXIT_SUCCESS;
}

// EOF
