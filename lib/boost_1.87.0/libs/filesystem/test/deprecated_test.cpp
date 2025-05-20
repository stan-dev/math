//  deprecated_test program --------------------------------------------------//

//  Copyright Beman Dawes 2002
//  Copyright Vladimir Prus 2002

//  Distributed under the Boost Software License, Version 1.0.
//  See http://www.boost.org/LICENSE_1_0.txt

//  Library home page: http://www.boost.org/libs/filesystem

//  This test verifies that various deprecated names still work. This is
//  important to preserve existing code that uses the old names.

#define BOOST_FILESYSTEM_DEPRECATED
#define BOOST_FILESYSTEM_ALLOW_DEPRECATED

#include <boost/filesystem.hpp>

#include <list>
#include <vector>

#include <boost/core/lightweight_test.hpp>
#include <boost/detail/lightweight_main.hpp>

namespace fs = boost::filesystem;
using boost::filesystem::path;

#define PATH_CHECK(a, b) check(a, b, __LINE__)

namespace {
std::string platform(BOOST_PLATFORM);

std::list< char > l;       // see main() for initialization to s, t, r, i, n, g
std::list< wchar_t > wl;   // see main() for initialization to w, s, t, r, i, n, g
std::vector< char > v;     // see main() for initialization to f, u, z
std::vector< wchar_t > wv; // see main() for initialization to w, f, u, z

void check(const fs::path& source, const std::string& expected, int line)
{
    if (source.generic_string() == expected)
        return;

    ++::boost::detail::test_errors();

    std::cout << __FILE__ << '(' << line << ") source.string(): \"" << source.string()
              << "\" != expected: \"" << expected
              << "\"" << std::endl;
}

//  path_container_ctor_test ---------------------------------------------------------//

void path_container_ctor_test()
{
    path x4v(v); // std::vector<char>
    PATH_CHECK(x4v, "fuz");
    BOOST_TEST_EQ(x4v.native().size(), 3U);

    path x5v(wv); // std::vector<wchar_t>
    PATH_CHECK(x5v, "wfuz");
    BOOST_TEST_EQ(x5v.native().size(), 4U);

    // non-contiguous containers
    path x10(l); // std::list<char>
    PATH_CHECK(x10, "string");
    BOOST_TEST_EQ(x10.native().size(), 6U);

    path xll(wl); // std::list<wchar_t>
    PATH_CHECK(xll, "wstring");
    BOOST_TEST_EQ(xll.native().size(), 7U);
}

} // unnamed namespace

//--------------------------------------------------------------------------------------//

int cpp_main(int /*argc*/, char* /*argv*/[])
{
    // The choice of platform is make at runtime rather than compile-time
    // so that compile errors for all platforms will be detected even though
    // only the current platform is runtime tested.
    platform = (platform == "Win32" || platform == "Win64" || platform == "Cygwin") ? "Windows" : "POSIX";
    std::cout << "Platform is " << platform << '\n';

    l.push_back('s');
    l.push_back('t');
    l.push_back('r');
    l.push_back('i');
    l.push_back('n');
    l.push_back('g');

    wl.push_back(L'w');
    wl.push_back(L's');
    wl.push_back(L't');
    wl.push_back(L'r');
    wl.push_back(L'i');
    wl.push_back(L'n');
    wl.push_back(L'g');

    v.push_back('f');
    v.push_back('u');
    v.push_back('z');

    wv.push_back(L'w');
    wv.push_back(L'f');
    wv.push_back(L'u');
    wv.push_back(L'z');

    BOOST_TEST(fs::initial_path() == fs::current_path());

    //path::default_name_check(fs::no_check);

    fs::directory_entry de("foo.bar", fs::file_status(fs::regular_file, fs::owner_all), fs::file_status(fs::directory_file, fs::group_all));

    BOOST_TEST(de.path() == "foo.bar");
    BOOST_TEST(de.status() == fs::file_status(fs::regular_file, fs::owner_all));
    BOOST_TEST(de.symlink_status() == fs::file_status(fs::directory_file, fs::group_all));
    BOOST_TEST(de < fs::directory_entry("goo.bar", fs::file_status(), fs::file_status()));
    BOOST_TEST(de == fs::directory_entry("foo.bar", fs::file_status(), fs::file_status()));
    BOOST_TEST(de != fs::directory_entry("goo.bar", fs::file_status(), fs::file_status()));
    de.replace_filename("bar.foo", fs::file_status(), fs::file_status());
    BOOST_TEST(de.path() == "bar.foo");

    const fs::path temp_dir(fs::current_path() / ".." / fs::unique_path("deprecated_test-%%%%-%%%%-%%%%"));
    std::cout << "temp_dir is " << temp_dir.string() << std::endl;

    fs::create_directory(temp_dir);

    path_container_ctor_test();

    std::cout << "post-test removal of " << temp_dir << std::endl;
    BOOST_TEST(fs::remove_all(temp_dir) != 0);

    return ::boost::report_errors();
}
