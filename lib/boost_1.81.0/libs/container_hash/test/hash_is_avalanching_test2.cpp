// Copyright 2022 Peter Dimov.
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/container_hash/hash.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <boost/config.hpp>
#include <boost/config/pragma_message.hpp>

#if defined(BOOST_NO_CXX17_HDR_STRING_VIEW)

BOOST_PRAGMA_MESSAGE( "Test skipped, BOOST_NO_CXX17_HDR_STRING_VIEW is defined" )
int main() {}

#else

#include <boost/unordered/hash_traits.hpp>
#include <string_view>

enum my_char { min = 0, max = 255 };

int main()
{
    using boost::unordered::hash_is_avalanching;

    BOOST_TEST_TRAIT_TRUE(( hash_is_avalanching< boost::hash<std::string_view> > ));
    BOOST_TEST_TRAIT_TRUE(( hash_is_avalanching< boost::hash<std::wstring_view> > ));

#if !defined(BOOST_NO_CXX11_CHAR16_T)

    BOOST_TEST_TRAIT_TRUE(( hash_is_avalanching< boost::hash<std::u16string_view> > ));

#endif

#if !defined(BOOST_NO_CXX11_CHAR32_T)

    BOOST_TEST_TRAIT_TRUE(( hash_is_avalanching< boost::hash<std::u32string_view> > ));

#endif

#if defined(__cpp_char8_t) && __cpp_char8_t >= 201811L

    BOOST_TEST_TRAIT_TRUE(( hash_is_avalanching< boost::hash< std::basic_string_view<char8_t> > > ));

#endif

    BOOST_TEST_TRAIT_FALSE(( hash_is_avalanching< boost::hash<std::basic_string_view<my_char> > > ));

    return boost::report_errors();
}

#endif
