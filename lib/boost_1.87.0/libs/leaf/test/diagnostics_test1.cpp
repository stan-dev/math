// Copyright 2018-2024 Emil Dotchevski and Reverge Studios, Inc.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef BOOST_LEAF_TEST_SINGLE_HEADER
#   include "leaf.hpp"
#else
#   include <boost/leaf/config.hpp>
#   include <boost/leaf/diagnostics.hpp>
#   include <boost/leaf/result.hpp>
#   include <boost/leaf/common.hpp>
#endif

#if BOOST_LEAF_CFG_STD_STRING
#   include <sstream>
#   include <iostream>
#endif

#include <algorithm>

#include "lightweight_test.hpp"

namespace leaf = boost::leaf;

enum class enum_class
{
    enum0,
    enum1
};

template <int A>
struct unexpected_test
{
    int value;
};

struct my_exception:
    virtual std::exception
{
    char const * what() const noexcept
    {
        return "my_exception what";
    }
};

struct printable_value
{
    friend std::ostream & operator<<( std::ostream & os, printable_value const & )
    {
        return os << "printed printable_value";
    }
};

struct non_printable_value
{
};

struct printable_info_printable_value
{
    printable_value value;

    friend std::ostream & operator<<( std::ostream & os, printable_info_printable_value const & x )
    {
        return os << x.value;
    }
};

struct printable_info_non_printable_value
{
    non_printable_value value;

    friend std::ostream & operator<<( std::ostream & os, printable_info_non_printable_value const & )
    {
        return os << "*** printable_info non_printable_value ***";
    }
};

class non_printable_info_printable_value
{
public:
    printable_value value;
};

struct non_printable_info_non_printable_value
{
    non_printable_value value;
};

struct hidden_printable_info_printable_value
{
    printable_value value;

    friend std::ostream & operator<<( std::ostream & os, hidden_printable_info_printable_value const & x )
    {
        return os << x.value << " ***";
    }
};

struct hidden_non_printable_info_printable_value
{
    printable_value value;
};

struct value_nullptr
{
    char const * value = nullptr;
};

struct value_ptr
{
    char const * value = "\"value\"";
};

namespace boost { namespace leaf {

template <> struct show_in_diagnostics<hidden_printable_info_printable_value>: std::false_type { };
template <> struct show_in_diagnostics<hidden_non_printable_info_printable_value>: std::false_type { };

} }

#if BOOST_LEAF_CFG_STD_STRING
namespace
{
    void remove_value(std::string & s, std::string const & search_prefix, std::string const & search_suffix)
    {
        auto search_prefix_size = search_prefix.size();
        auto s1 = s.find(search_prefix);
        if( s1 == s.npos )
            return;
        auto search_prefix_end = s1 + search_prefix_size;
        auto s2 = s.find(search_suffix, search_prefix_end);
        if( s2 == s.npos )
            return;
        s.replace(search_prefix_end, s2 - search_prefix_end, "<removed variance>");
    }

    void remove_variance(std::string & s)
    {
        s.erase(std::remove(s.begin(), s.end(), '\r'), s.end());
        remove_value(s, " reported at ", "\n");
        remove_value(s, "Caught:" BOOST_LEAF_CFG_DIAGNOSTICS_FIRST_DELIMITER, ": \"my_exception what\"");
    }

    bool cmp(std::string a, std::string const & b) noexcept
    {
        remove_variance(a);
        auto i = a.begin();
        auto j = b.begin();
        for( auto e = i + std::min(a.size(), b.size()); i != e; ++i, ++j )
            if( *i != *j )
            {
                std::cout <<
                    "a =\n----\n" << a << "----\n\n"
                    "b =\n----\n" << b << "----\n\n"
                    "Difference immediately after:\n" << std::string(a.begin(), i) << '\n';
                return false;
            }
        return i == a.end() && j == b.end();
    }
}
#endif

#define BOOST_LEAF_DIAGNOSTIC_INFO_NO_BOOST_LEAF_CFG_DIAGNOSTICS "\nboost::leaf::diagnostic_info N/A due to BOOST_LEAF_CFG_DIAGNOSTICS=0"
#define BOOST_LEAF_DIAGNOSTIC_DETAILS_NO_BOOST_LEAF_CFG_CAPTURE "\nboost::leaf::diagnostic_details N/A due to BOOST_LEAF_CFG_CAPTURE=0"
#define BOOST_LEAF_DIAGNOSTIC_DETAILS_NO_BOOST_LEAF_CFG_DIAGNOSTICS "\nboost::leaf::diagnostic_details N/A due to BOOST_LEAF_CFG_DIAGNOSTICS=0"

int main()
{
    std::cout << __LINE__  << " ---- result / error_info / no e_source_location\n";
    leaf::try_handle_all(
        []() -> leaf::result<void>
        {
            return BOOST_LEAF_NEW_ERROR(
                value_nullptr(),
                value_ptr(),
                42,
                printable_info_printable_value(),
                printable_info_non_printable_value(),
                non_printable_info_printable_value(),
                non_printable_info_non_printable_value(),
                hidden_printable_info_printable_value(),
                hidden_non_printable_info_printable_value(),
                enum_class::enum0,
                leaf::e_errno{ENOENT} );
        },
        [](
            value_nullptr,
            value_ptr,
            int,
            printable_info_printable_value,
            printable_info_non_printable_value,
            non_printable_info_printable_value,
            non_printable_info_non_printable_value,
            hidden_printable_info_printable_value,
            hidden_non_printable_info_printable_value,
            enum_class,
            leaf::e_errno,
            leaf::error_info const & unmatched )
        {
#if BOOST_LEAF_CFG_STD_STRING
            std::ostringstream st;
            st << unmatched;
            std::string s = st.str();
            std::cout << s << std::endl;
            BOOST_TEST(cmp(s,
                "Error with serial #1"
                "\n"
            ));

#endif
        },
        []()
        {
            BOOST_ERROR("Bad error dispatch");
        } );

    std::cout << __LINE__  << " ---- result / error_info\n";
    leaf::try_handle_all(
        []() -> leaf::result<void>
        {
            return BOOST_LEAF_NEW_ERROR(
                value_nullptr(),
                value_ptr(),
                42,
                printable_info_printable_value(),
                printable_info_non_printable_value(),
                non_printable_info_printable_value(),
                non_printable_info_non_printable_value(),
                hidden_printable_info_printable_value(),
                hidden_non_printable_info_printable_value(),
                enum_class::enum0,
                leaf::e_errno{ENOENT} );
        },
        [](
            value_nullptr,
            value_ptr,
            int,
            leaf::e_source_location,
            printable_info_printable_value,
            printable_info_non_printable_value,
            non_printable_info_printable_value,
            non_printable_info_non_printable_value,
            hidden_printable_info_printable_value,
            hidden_non_printable_info_printable_value,
            enum_class,
            leaf::e_errno,
            leaf::error_info const & unmatched )
        {
#if BOOST_LEAF_CFG_STD_STRING
            std::ostringstream st;
            st << unmatched;
            std::string s = st.str();
            std::cout << s << std::endl;
            BOOST_TEST(cmp(s,
                "Error with serial #2 reported at <removed variance>"
                "\n"
            ));

#endif
        },
        []()
        {
            BOOST_ERROR("Bad error dispatch");
        } );

    std::cout << __LINE__  << " ---- result / diagnostic_info\n";
    leaf::try_handle_all(
        []() -> leaf::result<void>
        {
            return BOOST_LEAF_NEW_ERROR(
                value_nullptr(),
                value_ptr(),
                42,
                printable_info_printable_value(),
                printable_info_non_printable_value(),
                non_printable_info_printable_value(),
                non_printable_info_non_printable_value(),
                hidden_printable_info_printable_value(),
                hidden_non_printable_info_printable_value(),
                enum_class::enum1,
                leaf::e_errno{ENOENT} );
        },
        [](
            value_nullptr,
            value_ptr,
            int,
            leaf::e_source_location,
            printable_info_printable_value,
            printable_info_non_printable_value,
            non_printable_info_printable_value,
            non_printable_info_non_printable_value,
            hidden_printable_info_printable_value,
            hidden_non_printable_info_printable_value,
            enum_class,
            leaf::e_errno,
            leaf::diagnostic_info const & unmatched )
        {
#if BOOST_LEAF_CFG_STD_STRING
            std::ostringstream st;
            st << unmatched;
            std::string s = st.str();
            std::cout << s << std::endl;
            if( BOOST_LEAF_CFG_DIAGNOSTICS )
                BOOST_TEST(cmp(s,
                    "Error with serial #3 reported at <removed variance>"
                    "\nCaught:"
                    BOOST_LEAF_CFG_DIAGNOSTICS_FIRST_DELIMITER "value_nullptr: <nullptr>"
                    BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "value_ptr: \"value\""
                    BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "int: 42"
                    BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "printable_info_printable_value: printed printable_value"
                    BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "printable_info_non_printable_value: *** printable_info non_printable_value ***"
                    BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "non_printable_info_printable_value: printed printable_value"
                    BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "non_printable_info_non_printable_value"
                    BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "enum_class: 1"
                    BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "boost::leaf::e_errno: 2, \"No such file or directory\""
                    "\n"
                ));
            else
                BOOST_TEST(cmp(s,
                    "Error with serial #3 reported at <removed variance>"
                    BOOST_LEAF_DIAGNOSTIC_INFO_NO_BOOST_LEAF_CFG_DIAGNOSTICS
                    "\n"
                ));
#endif
        },
        []()
        {
            BOOST_ERROR("Bad error dispatch");
        } );

    std::cout << __LINE__  << " ---- result / diagnostic_details\n";
    leaf::try_handle_all(
        []() -> leaf::result<void>
        {
            return BOOST_LEAF_NEW_ERROR(
                value_nullptr(),
                value_ptr(),
                42,
                printable_info_printable_value(),
                printable_info_non_printable_value(),
                non_printable_info_printable_value(),
                non_printable_info_non_printable_value(),
                hidden_printable_info_printable_value(),
                hidden_non_printable_info_printable_value(),
                unexpected_test<1>{1},
                unexpected_test<2>{2},
                enum_class::enum0,
                leaf::e_errno{ENOENT} );
        },
        [](
            value_nullptr,
            value_ptr,
            int,
            leaf::e_source_location,
            printable_info_printable_value,
            printable_info_non_printable_value,
            non_printable_info_printable_value,
            non_printable_info_non_printable_value,
            hidden_printable_info_printable_value,
            hidden_non_printable_info_printable_value,
            enum_class,
            leaf::e_errno,
            leaf::diagnostic_details const & di )
        {
#if BOOST_LEAF_CFG_STD_STRING
            std::ostringstream st;
            st << di;
            std::string s = st.str();
            std::cout << s << std::endl;
            if( BOOST_LEAF_CFG_DIAGNOSTICS )
                if( BOOST_LEAF_CFG_CAPTURE )
                    BOOST_TEST(cmp(s,
                        "Error with serial #4 reported at <removed variance>"
                        "\nCaught:"
                        BOOST_LEAF_CFG_DIAGNOSTICS_FIRST_DELIMITER "value_nullptr: <nullptr>"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "value_ptr: \"value\""
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "int: 42"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "printable_info_printable_value: printed printable_value"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "printable_info_non_printable_value: *** printable_info non_printable_value ***"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "non_printable_info_printable_value: printed printable_value"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "non_printable_info_non_printable_value"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "enum_class: 0"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "boost::leaf::e_errno: 2, \"No such file or directory\""
                        "\nDiagnostic details:"
                        BOOST_LEAF_CFG_DIAGNOSTICS_FIRST_DELIMITER "unexpected_test<1>: 1"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "unexpected_test<2>: 2"
                        "\n"
                    ));
                else
                    BOOST_TEST(cmp(s,
                        "Error with serial #4 reported at <removed variance>"
                        "\nCaught:"
                        BOOST_LEAF_CFG_DIAGNOSTICS_FIRST_DELIMITER "value_nullptr: <nullptr>"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "value_ptr: \"value\""
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "int: 42"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "printable_info_printable_value: printed printable_value"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "printable_info_non_printable_value: *** printable_info non_printable_value ***"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "non_printable_info_printable_value: printed printable_value"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "non_printable_info_non_printable_value"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "enum_class: 0"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "boost::leaf::e_errno: 2, \"No such file or directory\""
                        BOOST_LEAF_DIAGNOSTIC_DETAILS_NO_BOOST_LEAF_CFG_CAPTURE
                        "\n"
                    ));
            else
                BOOST_TEST(cmp(s,
                    "Error with serial #4 reported at <removed variance>"
                    BOOST_LEAF_DIAGNOSTIC_DETAILS_NO_BOOST_LEAF_CFG_DIAGNOSTICS
                    "\n"
                ));
#endif
        },
        []()
        {
            BOOST_ERROR("Bad error dispatch");
        } );

    std::cout << __LINE__  << " ---- result / diagnostic_details / nothing caught\n";
    leaf::try_handle_all(
        []() -> leaf::result<void>
        {
            return BOOST_LEAF_NEW_ERROR(
                value_nullptr(),
                value_ptr(),
                42,
                printable_info_printable_value(),
                printable_info_non_printable_value(),
                non_printable_info_printable_value(),
                non_printable_info_non_printable_value(),
                hidden_printable_info_printable_value(),
                hidden_non_printable_info_printable_value(),
                unexpected_test<1>{1},
                unexpected_test<2>{2},
                enum_class::enum0,
                leaf::e_errno{ENOENT} );
        },
        []( leaf::diagnostic_details const & di )
        {
#if BOOST_LEAF_CFG_STD_STRING
            std::ostringstream st;
            st << di;
            std::string s = st.str();
            std::cout << s << std::endl;
            if( BOOST_LEAF_CFG_DIAGNOSTICS )
                if( BOOST_LEAF_CFG_CAPTURE )
                    BOOST_TEST(cmp(s,
                        "Error with serial #5 reported at <removed variance>"
                        "\nDiagnostic details:"
                        BOOST_LEAF_CFG_DIAGNOSTICS_FIRST_DELIMITER "value_nullptr: <nullptr>"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "value_ptr: \"value\""
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "int: 42"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "printable_info_printable_value: printed printable_value"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "printable_info_non_printable_value: *** printable_info non_printable_value ***"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "non_printable_info_printable_value: printed printable_value"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "non_printable_info_non_printable_value"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "unexpected_test<1>: 1"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "unexpected_test<2>: 2"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "enum_class: 0"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "boost::leaf::e_errno: 2, \"No such file or directory\""
                        "\n"
                    ));
                else
                    BOOST_TEST(cmp(s,
                        "Error with serial #5 reported at <removed variance>"
                        BOOST_LEAF_DIAGNOSTIC_DETAILS_NO_BOOST_LEAF_CFG_CAPTURE
                        "\n"
                    ));
            else
                BOOST_TEST(cmp(s,
                    "Error with serial #5 reported at <removed variance>"
                    BOOST_LEAF_DIAGNOSTIC_DETAILS_NO_BOOST_LEAF_CFG_DIAGNOSTICS
                    "\n"
                ));
#endif
        },
        []()
        {
            BOOST_ERROR("Bad error dispatch");
        } );

    ////////////////////////////////////////

#ifndef BOOST_LEAF_NO_EXCEPTIONS

    std::cout << __LINE__  << " ---- exception / error_info / no e_source_location\n";
    leaf::try_catch(
        []
        {
            BOOST_LEAF_THROW_EXCEPTION( my_exception(),
                value_nullptr(),
                value_ptr(),
                42,
                printable_info_printable_value(),
                printable_info_non_printable_value(),
                non_printable_info_printable_value(),
                non_printable_info_non_printable_value(),
                hidden_printable_info_printable_value(),
                hidden_non_printable_info_printable_value(),
                enum_class::enum1,
                leaf::e_errno{ENOENT} );
        },
        [](
            value_nullptr,
            value_ptr,
            int,
            printable_info_printable_value,
            printable_info_non_printable_value,
            non_printable_info_printable_value,
            non_printable_info_non_printable_value,
            hidden_printable_info_printable_value,
            hidden_non_printable_info_printable_value,
            enum_class,
            leaf::e_errno,
            leaf::error_info const & unmatched )
        {
#if BOOST_LEAF_CFG_STD_STRING
            std::ostringstream st;
            st << unmatched;
            std::string s = st.str();
            std::cout << s << std::endl;
            BOOST_TEST(cmp(s,
                "Error with serial #6"
                "\nCaught:"
                BOOST_LEAF_CFG_DIAGNOSTICS_FIRST_DELIMITER "<removed variance>: \"my_exception what\""
                "\n"
            ));
#endif
        } );

    std::cout << __LINE__  << " ---- exception / error_info\n";
    leaf::try_catch(
        []
        {
            BOOST_LEAF_THROW_EXCEPTION( my_exception(),
                value_nullptr(),
                value_ptr(),
                42,
                printable_info_printable_value(),
                printable_info_non_printable_value(),
                non_printable_info_printable_value(),
                non_printable_info_non_printable_value(),
                hidden_printable_info_printable_value(),
                hidden_non_printable_info_printable_value(),
                enum_class::enum1,
                leaf::e_errno{ENOENT} );
        },
        [](
            value_nullptr,
            value_ptr,
            int,
            leaf::e_source_location,
            printable_info_printable_value,
            printable_info_non_printable_value,
            non_printable_info_printable_value,
            non_printable_info_non_printable_value,
            hidden_printable_info_printable_value,
            hidden_non_printable_info_printable_value,
            enum_class,
            leaf::e_errno,
            leaf::error_info const & unmatched )
        {
#if BOOST_LEAF_CFG_STD_STRING
            std::ostringstream st;
            st << unmatched;
            std::string s = st.str();
            std::cout << s << std::endl;
            BOOST_TEST(cmp(s,
                "Error with serial #7 reported at <removed variance>"
                "\nCaught:"
                BOOST_LEAF_CFG_DIAGNOSTICS_FIRST_DELIMITER "<removed variance>: \"my_exception what\""
                "\n"
            ));
#endif
        } );

    std::cout << __LINE__  << " ---- exception / diagnostic_info\n";
    leaf::try_catch(
        []
        {
            BOOST_LEAF_THROW_EXCEPTION( my_exception(),
                value_nullptr(),
                value_ptr(),
                42,
                printable_info_printable_value(),
                printable_info_non_printable_value(),
                non_printable_info_printable_value(),
                non_printable_info_non_printable_value(),
                hidden_printable_info_printable_value(),
                hidden_non_printable_info_printable_value(),
                enum_class::enum0,
                leaf::e_errno{ENOENT} );
        },
        [](
            value_nullptr,
            value_ptr,
            int,
            leaf::e_source_location,
            printable_info_printable_value,
            printable_info_non_printable_value,
            non_printable_info_printable_value,
            non_printable_info_non_printable_value,
            hidden_printable_info_printable_value,
            hidden_non_printable_info_printable_value,
            enum_class,
            leaf::e_errno,
            leaf::diagnostic_info const & unmatched )
        {
#if BOOST_LEAF_CFG_STD_STRING
            std::ostringstream st;
            st << unmatched;
            std::string s = st.str();
            std::cout << s << std::endl;
            if( BOOST_LEAF_CFG_DIAGNOSTICS )
                BOOST_TEST(cmp(s,
                    "Error with serial #8 reported at <removed variance>"
                    "\nCaught:"
                    BOOST_LEAF_CFG_DIAGNOSTICS_FIRST_DELIMITER "<removed variance>: \"my_exception what\""
                    BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "value_nullptr: <nullptr>"
                    BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "value_ptr: \"value\""
                    BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "int: 42"
                    BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "printable_info_printable_value: printed printable_value"
                    BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "printable_info_non_printable_value: *** printable_info non_printable_value ***"
                    BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "non_printable_info_printable_value: printed printable_value"
                    BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "non_printable_info_non_printable_value"
                    BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "enum_class: 0"
                    BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "boost::leaf::e_errno: 2, \"No such file or directory\""
                    "\n"
                ));
            else
                BOOST_TEST(cmp(s,
                    "Error with serial #8 reported at <removed variance>"
                    "\nCaught:"
                    BOOST_LEAF_CFG_DIAGNOSTICS_FIRST_DELIMITER "<removed variance>: \"my_exception what\""
                    BOOST_LEAF_DIAGNOSTIC_INFO_NO_BOOST_LEAF_CFG_DIAGNOSTICS
                    "\n"
                ));
#endif
        } );

    std::cout << __LINE__  << " ---- exception / diagnostic_details\n";
    leaf::try_catch(
        []
        {
            BOOST_LEAF_THROW_EXCEPTION( my_exception(),
                value_nullptr(),
                value_ptr(),
                42,
                printable_info_printable_value(),
                printable_info_non_printable_value(),
                non_printable_info_printable_value(),
                non_printable_info_non_printable_value(),
                hidden_printable_info_printable_value(),
                hidden_non_printable_info_printable_value(),
                enum_class::enum1,
                unexpected_test<1>{1},
                unexpected_test<2>{2},
                leaf::e_errno{ENOENT} );
        },
        [](
            value_nullptr,
            value_ptr,
            int,
            leaf::e_source_location,
            printable_info_printable_value,
            printable_info_non_printable_value,
            non_printable_info_printable_value,
            non_printable_info_non_printable_value,
            hidden_printable_info_printable_value,
            hidden_non_printable_info_printable_value,
            enum_class,
            leaf::e_errno,
            leaf::diagnostic_details const & di )
        {
#if BOOST_LEAF_CFG_STD_STRING
            std::ostringstream st;
            st << di;
            std::string s = st.str();
            std::cout << s << std::endl;
            if( BOOST_LEAF_CFG_DIAGNOSTICS )
                if( BOOST_LEAF_CFG_CAPTURE )
                    BOOST_TEST(cmp(s,
                        "Error with serial #9 reported at <removed variance>"
                        "\nCaught:"
                        BOOST_LEAF_CFG_DIAGNOSTICS_FIRST_DELIMITER "<removed variance>: \"my_exception what\""
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "value_nullptr: <nullptr>"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "value_ptr: \"value\""
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "int: 42"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "printable_info_printable_value: printed printable_value"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "printable_info_non_printable_value: *** printable_info non_printable_value ***"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "non_printable_info_printable_value: printed printable_value"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "non_printable_info_non_printable_value"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "enum_class: 1"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "boost::leaf::e_errno: 2, \"No such file or directory\""
                        "\nDiagnostic details:"
                        BOOST_LEAF_CFG_DIAGNOSTICS_FIRST_DELIMITER "unexpected_test<1>: 1"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "unexpected_test<2>: 2"
                        "\n"
                    ));
                else
                    BOOST_TEST(cmp(s,
                        "Error with serial #9 reported at <removed variance>"
                        "\nCaught:"
                        BOOST_LEAF_CFG_DIAGNOSTICS_FIRST_DELIMITER "<removed variance>: \"my_exception what\""
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "value_nullptr: <nullptr>"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "value_ptr: \"value\""
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "int: 42"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "printable_info_printable_value: printed printable_value"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "printable_info_non_printable_value: *** printable_info non_printable_value ***"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "non_printable_info_printable_value: printed printable_value"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "non_printable_info_non_printable_value"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "enum_class: 1"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "boost::leaf::e_errno: 2, \"No such file or directory\""
                        BOOST_LEAF_DIAGNOSTIC_DETAILS_NO_BOOST_LEAF_CFG_CAPTURE
                        "\n"
                    ));
            else
                BOOST_TEST(cmp(s,
                    "Error with serial #9 reported at <removed variance>"
                    "\nCaught:"
                    BOOST_LEAF_CFG_DIAGNOSTICS_FIRST_DELIMITER "<removed variance>: \"my_exception what\""
                    BOOST_LEAF_DIAGNOSTIC_DETAILS_NO_BOOST_LEAF_CFG_DIAGNOSTICS
                    "\n"
                ));
#endif
        } );

    std::cout << __LINE__  << " ---- exception / diagnostic_details / nothing caught\n";
    leaf::try_catch(
        []
        {
            BOOST_LEAF_THROW_EXCEPTION( my_exception(),
                value_nullptr(),
                value_ptr(),
                42,
                printable_info_printable_value(),
                printable_info_non_printable_value(),
                non_printable_info_printable_value(),
                non_printable_info_non_printable_value(),
                hidden_printable_info_printable_value(),
                hidden_non_printable_info_printable_value(),
                enum_class::enum1,
                unexpected_test<1>{1},
                unexpected_test<2>{2},
                leaf::e_errno{ENOENT} );
        },
        []( leaf::diagnostic_details const & di )
        {
#if BOOST_LEAF_CFG_STD_STRING
            std::ostringstream st;
            st << di;
            std::string s = st.str();
            std::cout << s << std::endl;
            if( BOOST_LEAF_CFG_DIAGNOSTICS )
                if( BOOST_LEAF_CFG_CAPTURE )
                    BOOST_TEST(cmp(s,
                        "Error with serial #10 reported at <removed variance>"
                        "\nCaught:"
                        BOOST_LEAF_CFG_DIAGNOSTICS_FIRST_DELIMITER "<removed variance>: \"my_exception what\""
                        "\nDiagnostic details:"
                        BOOST_LEAF_CFG_DIAGNOSTICS_FIRST_DELIMITER "value_nullptr: <nullptr>"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "value_ptr: \"value\""
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "int: 42"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "printable_info_printable_value: printed printable_value"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "printable_info_non_printable_value: *** printable_info non_printable_value ***"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "non_printable_info_printable_value: printed printable_value"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "non_printable_info_non_printable_value"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "enum_class: 1"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "unexpected_test<1>: 1"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "unexpected_test<2>: 2"
                        BOOST_LEAF_CFG_DIAGNOSTICS_DELIMITER "boost::leaf::e_errno: 2, \"No such file or directory\""
                        "\n"
                    ));
                else
                    BOOST_TEST(cmp(s,
                        "Error with serial #10 reported at <removed variance>"
                        "\nCaught:"
                        BOOST_LEAF_CFG_DIAGNOSTICS_FIRST_DELIMITER "<removed variance>: \"my_exception what\""
                        BOOST_LEAF_DIAGNOSTIC_DETAILS_NO_BOOST_LEAF_CFG_CAPTURE
                        "\n"
                    ));
            else
                BOOST_TEST(cmp(s,
                    "Error with serial #10 reported at <removed variance>"
                    "\nCaught:"
                    BOOST_LEAF_CFG_DIAGNOSTICS_FIRST_DELIMITER "<removed variance>: \"my_exception what\""
                    BOOST_LEAF_DIAGNOSTIC_DETAILS_NO_BOOST_LEAF_CFG_DIAGNOSTICS
                    "\n"
                ));
#endif
        } );

#endif

    return boost::report_errors();
}
