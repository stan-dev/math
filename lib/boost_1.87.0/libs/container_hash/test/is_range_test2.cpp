// Copyright 2022 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#if defined(_MSC_VER)
# pragma warning(disable: 4714) // forceinline not inlined
#endif

#if defined(__GNUC__) || defined(__clang__)
# pragma GCC diagnostic ignored "-Wsign-conversion"
#endif

#include <boost/container_hash/is_range.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <boost/config.hpp>
#include <boost/filesystem/path.hpp>

#if !defined(BOOST_NO_CXX17_HDR_FILESYSTEM) && !defined(__MINGW32__)
# include <filesystem>
#endif

int main()
{
    using boost::container_hash::is_range;

    BOOST_TEST_TRAIT_FALSE((is_range< boost::filesystem::path >));
    BOOST_TEST_TRAIT_FALSE((is_range< boost::filesystem::path const >));

#if !defined(BOOST_NO_CXX17_HDR_FILESYSTEM) && !defined(__MINGW32__)

    BOOST_TEST_TRAIT_FALSE((is_range< std::filesystem::path >));
    BOOST_TEST_TRAIT_FALSE((is_range< std::filesystem::path const >));

#endif

    return boost::report_errors();
}
