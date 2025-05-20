//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
// Copyright (c) 2022-2023 Alexander Grund
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/locale/encoding.hpp>
#include <boost/locale/generator.hpp>
#include <algorithm>

#include "boostLocale/test/tools.hpp"
#include "boostLocale/test/unit_test.hpp"

const bool test_iso_8859_8 =
#if defined(BOOST_LOCALE_WITH_ICU) || defined(BOOST_LOCALE_WITH_ICONV)
  true;
#else
  hasWinCodepage(28598);
#endif

#if defined(BOOST_LOCALE_WITH_ICONV)
// Reproduce issue #206 to detect faulty IConv
static bool isFaultyIconv()
{
    namespace blc = boost::locale::conv;
    auto from_utf = blc::detail::make_utf_decoder<char>("ISO-2022-CN", blc::skip, blc::detail::conv_backend::IConv);
    try {
        from_utf->convert("实");
    } catch(const std::runtime_error& e) {                                         // LCOV_EXCL_LINE
        return std::string(e.what()).find("IConv is faulty") != std::string::npos; // LCOV_EXCL_LINE
    }
    return false;
}
#else
constexpr bool isFaultyIconv()
{
    return false;
}
#endif

constexpr boost::locale::conv::detail::conv_backend all_conv_backends[] = {
#ifdef BOOST_LOCALE_WITH_ICONV
  boost::locale::conv::detail::conv_backend::IConv,
#endif
#ifdef BOOST_LOCALE_WITH_ICU
  boost::locale::conv::detail::conv_backend::ICU,
#endif
#if BOOST_LOCALE_USE_WIN32_API
  boost::locale::conv::detail::conv_backend::WinAPI,
#endif
};

std::ostream& operator<<(std::ostream& s, boost::locale::conv::detail::conv_backend impl)
{
    using boost::locale::conv::detail::conv_backend;
    switch(impl) {
        case conv_backend::Default: return s << "[Default]"; // LCOV_EXCL_LINE
        case conv_backend::IConv: return s << "[IConv]";
        case conv_backend::ICU: return s << "[ICU]";
        case conv_backend::WinAPI: return s << "[WinAPI]";
    }
    return s; // LCOV_EXCL_LINE
}

#define TEST_FAIL_CONVERSION(X) TEST_THROWS(X, boost::locale::conv::conversion_error)

template<typename Char>
void test_to_utf_for_impls(const std::string& source,
                           const std::basic_string<Char>& target,
                           const std::string& encoding,
                           const bool expectSuccess = true,
                           const bool test_default = true)
{
    if(test_default) {
        boost::locale::conv::utf_encoder<Char> conv(encoding);
        TEST_EQ(conv(source), target);
    }
    for(const auto impl : all_conv_backends) {
        std::cout << "----- " << impl << '\n';
        using boost::locale::conv::invalid_charset_error;
        try {
            auto convPtr =
              boost::locale::conv::detail::make_utf_encoder<Char>(encoding, boost::locale::conv::skip, impl);
            TEST_EQ(convPtr->convert(source), target);
        } catch(invalid_charset_error&) {
            std::cout << "--- Charset not supported\n"; // LCOV_EXCL_LINE
            continue;                                   // LCOV_EXCL_LINE
        }
        if(!expectSuccess) {
            auto convPtr =
              boost::locale::conv::detail::make_utf_encoder<Char>(encoding, boost::locale::conv::stop, impl);
            TEST_FAIL_CONVERSION(convPtr->convert(source));
        }
    }
    if(encoding == "UTF-8") {
        using boost::locale::conv::utf_to_utf;
        TEST_EQ(utf_to_utf<Char>(source), target);
        if(expectSuccess)
            TEST_EQ(utf_to_utf<char>(source), source);
        else
            TEST_FAIL_CONVERSION(utf_to_utf<Char>(source, boost::locale::conv::stop));
    }
}

template<typename Char>
void test_from_utf_for_impls(const std::basic_string<Char>& source,
                             const std::string& target,
                             const std::string& encoding,
                             const bool expectSuccess = true,
                             const bool test_default = true)
{
    if(test_default) {
        boost::locale::conv::utf_decoder<Char> conv(encoding);
        TEST_EQ(conv(source), target);
    }
    for(const auto impl : all_conv_backends) {
        std::cout << "----- " << impl << '\n';
        using boost::locale::conv::invalid_charset_error;
        try {
            auto convPtr =
              boost::locale::conv::detail::make_utf_decoder<Char>(encoding, boost::locale::conv::skip, impl);
            TEST_EQ(convPtr->convert(source), target);
        } catch(invalid_charset_error&) {
            std::cout << "--- Charset not supported\n"; // LCOV_EXCL_LINE
            continue;                                   // LCOV_EXCL_LINE
        }
        if(!expectSuccess) {
            auto convPtr =
              boost::locale::conv::detail::make_utf_decoder<Char>(encoding, boost::locale::conv::stop, impl);
            TEST_FAIL_CONVERSION(convPtr->convert(source));
        }
    }
    if(encoding == "UTF-8") {
        using boost::locale::conv::utf_to_utf;
        TEST_EQ(utf_to_utf<char>(source), target);
        if(expectSuccess)
            TEST_EQ(utf_to_utf<Char>(source), source);
        else
            TEST_FAIL_CONVERSION(utf_to_utf<char>(source, boost::locale::conv::stop));
    }
}

template<typename Char>
void test_to_from_utf(const std::string& source,
                      const std::basic_string<Char>& target,
                      const std::string& encoding,
                      const bool test_default = true)
{
    std::cout << "-- " << encoding << std::endl;

    if(test_default) {
        TEST_EQ(boost::locale::conv::to_utf<Char>(source, encoding), target);
        TEST_EQ(boost::locale::conv::from_utf<Char>(target, encoding), source);
    }
    test_to_utf_for_impls(source, target, encoding, true, test_default);
    test_from_utf_for_impls(target, source, encoding, true, test_default);
}

template<typename Char>
void test_error_to_utf(const std::string& source, const std::basic_string<Char>& target, const std::string& encoding)
{
    using boost::locale::conv::to_utf;
    using boost::locale::conv::stop;

    // Default: Replace, no error
    TEST_EQ(to_utf<Char>(source, encoding), target);
    // Test all overloads with method=stop -> error
    // source as string, C-String, range
    TEST_FAIL_CONVERSION(to_utf<Char>(source, encoding, stop));
    TEST_FAIL_CONVERSION(to_utf<Char>(source.c_str(), encoding, stop));
    TEST_FAIL_CONVERSION(to_utf<Char>(source.c_str(), source.c_str() + source.size(), encoding, stop));
    // Same but encoding via locale
    const std::locale l = boost::locale::generator{}("en_US." + encoding);
    TEST_FAIL_CONVERSION(to_utf<Char>(source, l, stop));
    TEST_FAIL_CONVERSION(to_utf<Char>(source.c_str(), l, stop));
    TEST_FAIL_CONVERSION(to_utf<Char>(source.c_str(), source.c_str() + source.size(), l, stop));
    test_to_utf_for_impls(source, target, encoding, false);
}

template<typename Char>
void test_error_from_utf(const std::basic_string<Char>& source, const std::string& target, const std::string& encoding)
{
    using boost::locale::conv::from_utf;
    using boost::locale::conv::stop;

    // Default: Replace, no error
    TEST_EQ(from_utf<Char>(source, encoding), target);
    // Test all overloads with method=stop -> error
    // source as string, C-String, range
    TEST_FAIL_CONVERSION(from_utf<Char>(source, encoding, stop));
    TEST_FAIL_CONVERSION(from_utf<Char>(source.c_str(), encoding, stop));
    TEST_FAIL_CONVERSION(from_utf<Char>(source.c_str(), source.c_str() + source.size(), encoding, stop));
    // Same but encoding via locale
    const std::locale l = boost::locale::generator{}("en_US." + encoding);
    TEST_FAIL_CONVERSION(from_utf<Char>(source, l, stop));
    TEST_FAIL_CONVERSION(from_utf<Char>(source.c_str(), l, stop));
    TEST_FAIL_CONVERSION(from_utf<Char>(source.c_str(), source.c_str() + source.size(), l, stop));
    test_from_utf_for_impls(source, target, encoding, false);
}

template<typename Char>
std::basic_string<Char> utf(const std::string& s)
{
    return to<Char>(s);
}

template<>
std::basic_string<char> utf(const std::string& s)
{
    return s;
}

template<typename Char>
void test_with_0()
{
    std::cout << "-- Test string containing NULL chars" << std::endl;
    const char with_null[] = "foo\0\0 of\0";
    const std::string s_with_null(with_null, sizeof(with_null) - 1);
    const std::basic_string<Char> s_with_null2 = ascii_to<Char>(with_null);
    for(const std::string charset : {"UTF-8", "ISO8859-1"}) {
        for(const auto impl : all_conv_backends) {
            std::cout << "--- " << charset << " to UTF with Impl " << impl << std::endl;
            auto to_utf =
              boost::locale::conv::detail::make_utf_encoder<Char>(charset, boost::locale::conv::default_method, impl);
            TEST_EQ(to_utf->convert(s_with_null), s_with_null2);
            std::cout << "--- " << charset << " from UTF with Impl " << impl << std::endl;
            auto from_utf =
              boost::locale::conv::detail::make_utf_decoder<Char>(charset, boost::locale::conv::default_method, impl);
            TEST_EQ(from_utf->convert(s_with_null2), s_with_null);
        }
    }
    using boost::locale::conv::utf_to_utf;
    TEST_EQ(utf_to_utf<Char>(s_with_null), s_with_null2);
    TEST_EQ(utf_to_utf<Char>(s_with_null2), s_with_null2);
    TEST_EQ(utf_to_utf<char>(s_with_null2), s_with_null);
    TEST_EQ(utf_to_utf<char>(s_with_null), s_with_null);
}

template<typename Char, int n = sizeof(Char)>
struct utfutf;

#ifdef BOOST_MSVC
#    pragma warning(push)
#    pragma warning(disable : 4309) // narrowing static_cast warning
#endif
template<typename U8Char>
struct utfutf<U8Char, 1> {
    static const U8Char* ok() { return reinterpret_cast<const U8Char*>("grüßen"); }
    static const U8Char* bad()
    {
        return reinterpret_cast<const U8Char*>("gr\xFF"
                                               "üßen");
        // split into 2 to make SunCC happy
    }
    static U8Char bad_char() { return static_cast<U8Char>(0xFF); }
};

template<>
struct utfutf<wchar_t, 2> {
    static const wchar_t* ok() { return L"\x67\x72\xfc\xdf\x65\x6e"; }
    static const wchar_t* bad()
    {
        static wchar_t buf[256] = L"\x67\x72\xFF\xfc\xFE\xFD\xdf\x65\x6e";
        buf[2] = 0xDC01; // second surrogate must not be
        buf[4] = 0xD801; // First
        buf[5] = 0xD801; // Must be surrogate trail
        return buf;
    }
    static wchar_t bad_char() { return static_cast<wchar_t>(0xDC01); }
};

template<>
struct utfutf<wchar_t, 4> {
    static const wchar_t* ok() { return L"\x67\x72\xfc\xdf\x65\x6e"; }
    static const wchar_t* bad()
    {
        static wchar_t buf[256] = L"\x67\x72\xFF\xfc\xdf\x65\x6e";
        buf[2] = static_cast<wchar_t>(0x1000000); // > 10FFFF
        return buf;
    }
    static wchar_t bad_char() { return static_cast<wchar_t>(0x1000000); }
};
#ifdef BOOST_MSVC
#    pragma warning(pop)
#endif

template<typename CharOut, typename CharIn>
void test_combinations()
{
    using boost::locale::conv::utf_to_utf;
    typedef utfutf<CharOut> out;
    typedef utfutf<CharIn> in;
    const CharIn* inOk = in::ok();
    // Both overloads: C-string and string. Both call the range overload
    TEST((utf_to_utf<CharOut>(inOk) == out::ok()));
    TEST((utf_to_utf<CharOut>(std::basic_string<CharIn>(inOk)) == out::ok()));
    const CharIn* inBad = in::bad();
    // Again both overloads
    TEST_FAIL_CONVERSION((utf_to_utf<CharOut>(inBad, boost::locale::conv::stop)));
    TEST_FAIL_CONVERSION((utf_to_utf<CharOut>(std::basic_string<CharIn>(inBad), boost::locale::conv::stop)));
    TEST((utf_to_utf<CharOut>(in::bad()) == out::ok()));
}

void test_all_combinations()
{
    std::cout << "Testing utf_to_utf\n";
    std::cout << "  char<-char" << std::endl;
    test_combinations<char, char>();
    std::cout << "  char<-wchar" << std::endl;
    test_combinations<char, wchar_t>();
    std::cout << "  wchar<-char" << std::endl;
    test_combinations<wchar_t, char>();
    std::cout << "  wchar<-wchar" << std::endl;
    test_combinations<wchar_t, wchar_t>();
}

template<typename Char>
void test_utf_for()
{
    using boost::locale::conv::invalid_charset_error;

    {
        using boost::locale::conv::to_utf;
        using boost::locale::conv::from_utf;
        TEST_THROWS(to_utf<Char>("Hello", "invalid-charset"), invalid_charset_error);
        TEST_THROWS(from_utf<Char>(ascii_to<Char>("Hello"), "invalid-charset"), invalid_charset_error);
    }

    test_to_from_utf<Char>(to<char>("grüßen"), utf<Char>("grüßen"), "ISO8859-1");
    if(test_iso_8859_8)
        test_to_from_utf<Char>("\xf9\xec\xe5\xed", utf<Char>("שלום"), "ISO8859-8");
    test_to_from_utf<Char>("grüßen", utf<Char>("grüßen"), "UTF-8");
    test_to_from_utf<Char>("abc\"\xf0\xa0\x82\x8a\"", utf<Char>("abc\"\xf0\xa0\x82\x8a\""), "UTF-8");
    // Testing a codepage which may be an issue on Windows, see issue #121
    try {
        test_to_from_utf<Char>("\x1b$BE_5(\x1b(B", utf<Char>("冬季"), "iso-2022-jp");
    } catch(const invalid_charset_error&) { // LCOV_EXCL_LINE
        std::cout << "--- not supported\n"; // LCOV_EXCL_LINE
    }
    if(!isFaultyIconv()) {
        // Testing a codepage which may crash with IConv on macOS, see issue #196
        test_to_from_utf<Char>("\xa1\xad\xa1\xad", utf<Char>("……"), "gbk", false);
        // This might cause a bogus E2BIG on macOS, see issue #206
        test_to_from_utf<Char>("\x1b\x24\x29\x41\x0e\x4a\x35\xf", utf<Char>("实"), "ISO-2022-CN", false);
    }

    std::cout << "- Testing correct invalid bytes skipping\n";
    {
        std::cout << "-- UTF-8" << std::endl;

        std::cout << "--- At start single" << std::endl;
        test_error_to_utf<Char>("\xFFgrüßen", utf<Char>("grüßen"), "UTF-8");
        std::cout << "--- At start multiple" << std::endl;
        test_error_to_utf<Char>("\xFF\xFFgrüßen", utf<Char>("grüßen"), "UTF-8");

        std::cout << "--- At middle single" << std::endl;
        test_error_to_utf<Char>("g\xFFrüßen", utf<Char>("grüßen"), "UTF-8");
        std::cout << "--- At middle multiple" << std::endl;
        test_error_to_utf<Char>("g\xFF\xFF\xFFrüßen", utf<Char>("grüßen"), "UTF-8");

        std::cout << "--- At end single" << std::endl;
        test_error_to_utf<Char>("grüßen\xFF", utf<Char>("grüßen"), "UTF-8");
        std::cout << "--- At end multiple" << std::endl;
        test_error_to_utf<Char>("grüßen\xFF\xFF", utf<Char>("grüßen"), "UTF-8");

        try {
            std::cout << "-- ISO-8859-8" << std::endl;
            test_error_to_utf<Char>("\xFB", utf<Char>(""), "ISO-8859-8");
            test_error_to_utf<Char>("\xFB-", utf<Char>("-"), "ISO-8859-8");
            test_error_to_utf<Char>("test \xE0\xE1\xFB", utf<Char>("test \xd7\x90\xd7\x91"), "ISO-8859-8");
            test_error_to_utf<Char>("test \xE0\xE1\xFB-", utf<Char>("test \xd7\x90\xd7\x91-"), "ISO-8859-8");
        } catch(const invalid_charset_error&) { // LCOV_EXCL_LINE
            std::cout << "--- not supported\n"; // LCOV_EXCL_LINE
        }
        try {
            std::cout << "-- cp932" << std::endl;
            test_error_to_utf<Char>("\x83\xF8", utf<Char>(""), "cp932");
            test_error_to_utf<Char>("\x83\xF8-", utf<Char>("-"), "cp932");
            test_error_to_utf<Char>("test\xE0\xA0 \x83\xF8", utf<Char>("test\xe7\x87\xbf "), "cp932");
            test_error_to_utf<Char>("test\xE0\xA0 \x83\xF8-", utf<Char>("test\xe7\x87\xbf -"), "cp932");
        } catch(const invalid_charset_error&) { // LCOV_EXCL_LINE
            std::cout << "--- not supported\n"; // LCOV_EXCL_LINE
        }
        std::cout << "-- Error for encoding at start" << std::endl;
        test_error_from_utf<Char>(utf<Char>("שלום hello"), " hello", "ISO8859-1");
        std::cout << "-- Error for encoding at  middle and end" << std::endl;
        test_error_from_utf<Char>(utf<Char>("hello שלום world"), "hello  world", "ISO8859-1");
        std::cout << "-- Error for encoding at end" << std::endl;
        test_error_from_utf<Char>(utf<Char>("hello שלום"), "hello ", "ISO8859-1");
        std::cout << "-- Error for decoding to UTF-8" << std::endl;
        test_error_from_utf<Char>(utfutf<Char>::bad(), utfutf<char>::ok(), "UTF-8");
        std::cout << "-- Error for decoding to Latin1" << std::endl;
        test_error_from_utf<Char>(utfutf<Char>::bad(), to<char>(utfutf<char>::ok()), "Latin1");

        const std::basic_string<Char> onlyInvalidUtf(2, utfutf<Char>::bad_char());
        std::cout << "-- Error decoding string of only invalid chars to UTF-8" << std::endl;
        test_error_from_utf<Char>(onlyInvalidUtf, "", "UTF-8");
        std::cout << "-- Error decoding string of only invalid chars to Latin1" << std::endl;
        test_error_from_utf<Char>(onlyInvalidUtf, "", "Latin1");
    }

    test_with_0<Char>();
}

template<typename Char1, typename Char2>
void test_utf_to_utf_for(const std::string& utf8_string)
{
    const auto utf_string1 = utf<Char1>(utf8_string);
    const auto utf_string2 = utf<Char2>(utf8_string);
    using boost::locale::conv::utf_to_utf;
    TEST_EQ(utf_to_utf<Char1>(utf_string2), utf_string1);
    TEST_EQ(utf_to_utf<Char2>(utf_string1), utf_string2);
    TEST_EQ(utf_to_utf<Char1>(utf_string1), utf_string1);
    TEST_EQ(utf_to_utf<Char2>(utf_string2), utf_string2);
}

template<typename Char>
void test_utf_to_utf_for()
{
    const std::string& utf8_string = "A-Za-z0-9grüße'\xf0\xa0\x82\x8a'\xf4\x8f\xbf\xbf";
    std::cout << "---- char\n";
    test_utf_to_utf_for<Char, char>(utf8_string);
    test_to_utf_for_impls(utf8_string, utf<Char>(utf8_string), "UTF-8");
    test_from_utf_for_impls(utf<Char>(utf8_string), utf8_string, "UTF-8");
    std::cout << "---- wchar_t\n";
    test_utf_to_utf_for<Char, wchar_t>(utf8_string);
#ifndef BOOST_LOCALE_NO_CXX20_STRING8
    std::cout << "---- char8_t\n";
    test_utf_to_utf_for<Char, char8_t>(utf8_string);
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR16_T
    std::cout << "---- char16_t\n";
    test_utf_to_utf_for<Char, char16_t>(utf8_string);
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR32_T
    std::cout << "---- char32_t\n";
    test_utf_to_utf_for<Char, char32_t>(utf8_string);
#endif
}

void test_utf_to_utf()
{
    std::cout << "- Testing UTF to UTF conversion\n";
    std::cout << "-- char\n";
    test_utf_to_utf_for<char>();
    std::cout << "-- wchar_t\n";
    test_utf_to_utf_for<wchar_t>();
#ifndef BOOST_LOCALE_NO_CXX20_STRING8
    std::cout << "-- char8_t\n";
    test_utf_to_utf_for<char8_t>();
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR16_T
    std::cout << "-- char16_t\n";
    test_utf_to_utf_for<char16_t>();
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR32_T
    std::cout << "-- char32_t\n";
    test_utf_to_utf_for<char32_t>();
#endif
}

/// Allocator that reports when it has been used in a static variable
int globalUsedId = 0;
template<typename T>
struct CustomAllocator {
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using reference = T&;
    using const_reference = const T&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using propagate_on_container_move_assignment = std::true_type;
    using is_always_equal = std::false_type;

    template<typename U>
    struct rebind {
        typedef CustomAllocator<U> other;
    };

    CustomAllocator(const int id = 1) : id(id) {}
    template<typename U>
    CustomAllocator(const CustomAllocator<U>& other) : id(other.id)
    {}

    T* allocate(size_t n)
    {
        // Only count allocations of (w)chars, not e.g. internal proxy instances
        BOOST_LOCALE_START_CONST_CONDITION
        if(std::is_same<T, char>::value || std::is_same<T, wchar_t>::value)
            usedId += id;
        BOOST_LOCALE_END_CONST_CONDITION
        return base.allocate(n);
    }

    void deallocate(T* p, size_t n) { return base.deallocate(p, n); }

    static int& usedId;
    int id;

private:
    std::allocator<T> base;
};
template<class T, class U>
bool operator==(const CustomAllocator<T>&, const CustomAllocator<U>&)
{
    return true;
}

template<class T, class U>
bool operator!=(const CustomAllocator<T>&, const CustomAllocator<U>&)
{
    return false;
}

namespace detail {
// Note that using a static class variable does not work due to possible rebinds
int allocUsedId = 0;
} // namespace detail
template<typename T>
int& CustomAllocator<T>::usedId = detail::allocUsedId;

void test_utf_to_utf_allocator_support()
{
    using Alloc = CustomAllocator<wchar_t>;
    using AllocIn = CustomAllocator<char>;
    using boost::locale::conv::utf_to_utf;
    const auto method = boost::locale::conv::default_method;
    const std::string input(65, '0'); // Long enough to avoid SBO
    const AllocIn inputAllocator(17);
    const std::basic_string<char, std::char_traits<char>, AllocIn> inputWithAlloc(input.begin(),
                                                                                  input.end(),
                                                                                  inputAllocator);
    const std::basic_string<wchar_t, std::char_traits<wchar_t>, Alloc> output(input.begin(), input.end());
    const char* sBegin = input.data();
    const char* sEnd = sBegin + input.size();

    // Allocator via template param
    Alloc::usedId = 0;
    TEST_EQ((utf_to_utf<wchar_t, char, Alloc>(sBegin, sEnd)), output);
    TEST_EQ(Alloc::usedId, 1);
    Alloc::usedId = 0;
    TEST_EQ((utf_to_utf<wchar_t, char, Alloc>(sBegin, method)), output);
    TEST_EQ(Alloc::usedId, 1);
    Alloc::usedId = 0;
    TEST_EQ((utf_to_utf<wchar_t, char, Alloc>(inputWithAlloc)), output);
    TEST_EQ(Alloc::usedId, 1);
    Alloc::usedId = 0;
    TEST_EQ((utf_to_utf<wchar_t, char, Alloc>(inputWithAlloc, method)), output);
    TEST_EQ(Alloc::usedId, 1);

    // Pass allocator explicitly
    Alloc::usedId = 0;
    TEST_EQ(utf_to_utf<wchar_t>(sBegin, sEnd, method, Alloc(2)), output);
    TEST_EQ(Alloc::usedId, 2);
    Alloc::usedId = 0;
    TEST_EQ(utf_to_utf<wchar_t>(sBegin, method, Alloc(3)), output);
    TEST_EQ(Alloc::usedId, 3);
    Alloc::usedId = 0;
    TEST_EQ(utf_to_utf<wchar_t>(inputWithAlloc, method, Alloc(4)), output);
    TEST_EQ(Alloc::usedId, 4);
    // Same with using the default method
    Alloc::usedId = 0;
    TEST_EQ(utf_to_utf<wchar_t>(sBegin, sEnd, Alloc(2)), output);
    TEST_EQ(Alloc::usedId, 2);
    Alloc::usedId = 0;
    TEST_EQ(utf_to_utf<wchar_t>(sBegin, Alloc(3)), output);
    TEST_EQ(Alloc::usedId, 3);
    Alloc::usedId = 0;
    TEST_EQ(utf_to_utf<wchar_t>(inputWithAlloc, Alloc(4)), output);
    TEST_EQ(Alloc::usedId, 4);

    // Use allocator from input
    Alloc::usedId = 0;
    TEST_EQ(utf_to_utf<wchar_t>(inputWithAlloc), output);
    TEST_EQ(Alloc::usedId, inputAllocator.id);
    Alloc::usedId = 0;
    TEST_EQ(utf_to_utf<wchar_t>(inputWithAlloc, method), output);
    TEST_EQ(Alloc::usedId, inputAllocator.id);

    // Unchanged allocator for string overloads to check for ambiguous overloads
    AllocIn::usedId = 0;
    TEST_EQ(utf_to_utf<char>(inputWithAlloc, method, AllocIn(4)), inputWithAlloc);
    TEST_EQ(AllocIn::usedId, 4);
    TEST_EQ(utf_to_utf<char>(inputWithAlloc), inputWithAlloc);
    TEST_EQ(AllocIn::usedId, 4 + inputAllocator.id);
}

/// Test all overloads of to_utf/from_utf templated by Char
template<typename Char>
void test_latin1_conversions_for()
{
    const std::string utf8_string = "A-Za-z0-9grüße";
    const std::string sLatin1 = to<char>(utf8_string);
    // Sanity check that utf8_string is UTF-8 encoded (using multiple bytes for the special chars)
    // and sLatin1 is not encoded (1 byte per char)
    TEST_GT(utf8_string.length(), sLatin1.length());
    const std::basic_string<Char> sWide = utf<Char>(utf8_string);
    const std::string encoding = "Latin1";

    using boost::locale::conv::to_utf;
    using boost::locale::conv::utf_encoder;
    // 3 variants for source: string, C-string, range
    TEST_EQ(to_utf<Char>(sLatin1, encoding), sWide);
    TEST_EQ(to_utf<Char>(sLatin1.c_str(), encoding), sWide);
    TEST_EQ(to_utf<Char>(sLatin1.c_str(), sLatin1.c_str() + sLatin1.size(), encoding), sWide);
    TEST_EQ(utf_encoder<Char>(encoding)(sLatin1), sWide);
    TEST_EQ(utf_encoder<Char>(encoding).convert(sLatin1), sWide);
    TEST_EQ(utf_encoder<Char>(encoding).convert(sLatin1.c_str(), sLatin1.c_str() + sLatin1.size()), sWide);
    // Same but encoding given via locale
    const std::locale l = boost::locale::generator{}("en_US.Latin1");
    TEST_EQ(to_utf<Char>(sLatin1, l), sWide);
    TEST_EQ(to_utf<Char>(sLatin1.c_str(), l), sWide);
    TEST_EQ(to_utf<Char>(sLatin1.c_str(), sLatin1.c_str() + sLatin1.size(), l), sWide);

    using boost::locale::conv::from_utf;
    using boost::locale::conv::utf_decoder;
    // 3 variants for source: string, C-string, range
    TEST_EQ(from_utf<Char>(sWide, encoding), sLatin1);
    TEST_EQ(from_utf<Char>(sWide.c_str(), encoding), sLatin1);
    TEST_EQ(from_utf<Char>(sWide.c_str(), sWide.c_str() + sWide.size(), encoding), sLatin1);
    TEST_EQ(utf_decoder<Char>(encoding)(sWide), sLatin1);
    TEST_EQ(utf_decoder<Char>(encoding).convert(sWide), sLatin1);
    TEST_EQ(utf_decoder<Char>(encoding).convert(sWide.c_str(), sWide.c_str() + sWide.size()), sLatin1);
    // Same but encoding given via locale
    TEST_EQ(from_utf<Char>(sWide, l), sLatin1);
    TEST_EQ(from_utf<Char>(sWide.c_str(), l), sLatin1);
    TEST_EQ(from_utf<Char>(sWide.c_str(), sWide.c_str() + sWide.size(), l), sLatin1);

    // Empty string doesn't error/assert
    TEST_EQ(to_utf<Char>("", encoding), utf<Char>(""));
    TEST_EQ(from_utf<Char>(utf<Char>(""), encoding), std::string());
    test_to_utf_for_impls("", utf<Char>(""), encoding);
    test_from_utf_for_impls(utf<Char>(""), "", encoding);
}

/// Quick check of to_utf/from_utf overloads using the simple Latin1 encoding
void test_latin1_conversions()
{
    std::cout << "- Testing Latin1 conversion\n";
    std::cout << "-- char\n";
    test_latin1_conversions_for<char>();
    std::cout << "-- wchar_t\n";
    test_latin1_conversions_for<wchar_t>();
#ifndef BOOST_LOCALE_NO_CXX20_STRING8
    std::cout << "-- char8_t\n";
    test_latin1_conversions_for<char8_t>();
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR16_T
    std::cout << "-- char16_t\n";
    test_latin1_conversions_for<char16_t>();
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR32_T
    std::cout << "-- char32_t\n";
    test_latin1_conversions_for<char32_t>();
#endif
}

void test_between_for_impls(const std::string& source,
                            const std::string& target,
                            const std::string& to_encoding,
                            const std::string& from_encoding,
                            const bool expectSuccess = true)
{
    boost::locale::conv::narrow_converter conv(from_encoding, to_encoding);
    TEST_EQ(conv(source), target);
    for(const auto impl : all_conv_backends) {
        using boost::locale::conv::detail::make_narrow_converter;
        std::cout << "----- " << impl << '\n';
        using boost::locale::conv::invalid_charset_error;
        try {
            auto convPtr = make_narrow_converter(from_encoding, to_encoding, boost::locale::conv::skip, impl);
            TEST_EQ(convPtr->convert(source), target);
        } catch(invalid_charset_error&) {
            continue; // LCOV_EXCL_LINE
        }
        if(!expectSuccess) {
            auto convPtr = make_narrow_converter(from_encoding, to_encoding, boost::locale::conv::stop, impl);
            TEST_FAIL_CONVERSION(convPtr->convert(source));
        }
    }
    if(to_encoding == "UTF-8" && from_encoding == "UTF-8") {
        using boost::locale::conv::utf_to_utf;
        TEST_EQ(utf_to_utf<char>(source), target);
        if(!expectSuccess)
            TEST_FAIL_CONVERSION(utf_to_utf<char>(source, boost::locale::conv::stop));
    }
}

void test_error_between(const std::string& source,
                        const std::string& target,
                        const std::string& to_encoding,
                        const std::string& from_encoding)
{
    using boost::locale::conv::between;
    TEST_EQ(between(source, to_encoding, from_encoding), target);
    using boost::locale::conv::stop;
    TEST_FAIL_CONVERSION(between(source, to_encoding, from_encoding, stop));
    TEST_FAIL_CONVERSION(between(source.c_str(), to_encoding, from_encoding, stop));
    TEST_FAIL_CONVERSION(between(source.c_str(), source.c_str() + source.size(), to_encoding, from_encoding, stop));
    test_between_for_impls(source, target, to_encoding, from_encoding, false);
}

void test_between()
{
    using boost::locale::conv::between;
    const std::string utf8_string = "A-Za-z0-9grüße";
    const std::string sLatin1 = to<char>(utf8_string);
    TEST_GT(utf8_string.length(), sLatin1.length()); // Assert UTF encoding -> multi byte
    TEST_EQ(between(sLatin1, "UTF-8", "Latin1"), utf8_string);
    TEST_EQ(between(sLatin1.c_str(), "UTF-8", "Latin1"), utf8_string);
    TEST_EQ(between(sLatin1.c_str(), sLatin1.c_str() + sLatin1.size(), "UTF-8", "Latin1"), utf8_string);
    test_between_for_impls(sLatin1, utf8_string, "UTF-8", "Latin1");
    TEST_EQ(between(utf8_string, "Latin1", "UTF-8"), sLatin1);
    TEST_EQ(between(utf8_string.c_str(), "Latin1", "UTF-8"), sLatin1);
    TEST_EQ(between(utf8_string.c_str(), utf8_string.c_str() + utf8_string.size(), "Latin1", "UTF-8"), sLatin1);
    test_between_for_impls(utf8_string, sLatin1, "Latin1", "UTF-8");
    // Same encoding
    TEST_EQ(between(utf8_string, "UTF-8", "UTF-8"), utf8_string);
    test_between_for_impls(utf8_string, utf8_string, "UTF-8", "UTF-8");
    TEST_EQ(between(sLatin1, "Latin1", "Latin1"), sLatin1);
    test_between_for_impls(sLatin1, sLatin1, "Latin1", "Latin1");
    // Wrong encoding throws
    {
        using boost::locale::conv::invalid_charset_error;
        TEST_THROWS(between(sLatin1, "Invalid-Encoding", "Latin1"), invalid_charset_error);
        TEST_THROWS(between(sLatin1, "UTF-8", "Invalid-Encoding"), invalid_charset_error);
        TEST_THROWS(between(sLatin1, "Invalid-Encoding", "Invalid-Encoding"), invalid_charset_error);
        for(const auto impl : all_conv_backends) {
            std::cout << "----- " << impl << '\n';
            using boost::locale::conv::invalid_charset_error;
            using boost::locale::conv::skip;
            using boost::locale::conv::detail::make_narrow_converter;
            TEST_THROWS(make_narrow_converter("Invalid-Encoding", "Latin1", skip, impl), invalid_charset_error);
            TEST_THROWS(make_narrow_converter("UTF-8", "Invalid-Encoding", skip, impl), invalid_charset_error);
            TEST_THROWS(make_narrow_converter("Invalid-Encoding", "Invalid-Encoding", skip, impl),
                        invalid_charset_error);
        }
    }
    // Error handling
    // Unencodable char at start, middle, end
    test_error_between("שלום hello", " hello", "ISO8859-1", "UTF-8");
    test_error_between("hello שלום world", "hello  world", "ISO8859-1", "UTF-8");
    test_error_between("hello שלום", "hello ", "ISO8859-1", "UTF-8");
    // Undecodable char(s) at start, middle, end
    test_error_between("\xFFxfoo", "xfoo", "ISO8859-1", "UTF-8");
    test_error_between("\xFF\xFFyfoo", "yfoo", "ISO8859-1", "UTF-8");
    test_error_between("f\xFFoo2", "foo2", "ISO8859-1", "UTF-8");
    test_error_between("f\xFF\xFF\xFFoo3", "foo3", "ISO8859-1", "UTF-8");
    test_error_between("foo4\xFF", "foo4", "ISO8859-1", "UTF-8");
    test_error_between("foo5\xFF\xFF", "foo5", "ISO8859-1", "UTF-8");
    // Same but UTF-8 to UTF-8
    test_error_between("\xFFzfoo", "zfoo", "UTF-8", "UTF-8");
    test_error_between("f\xFFoo6", "foo6", "UTF-8", "UTF-8");
    test_error_between("f\xFF\xFF\xFFoo7", "foo7", "UTF-8", "UTF-8");
}

void test_utf_name();
void test_simple_encodings();
void test_win_codepages();

void test_main(int /*argc*/, char** /*argv*/)
{
    // Sanity check to<char>
    TEST_EQ(to<char>("grüßen"),
            "gr\xFC\xDF"
            "en");
    TEST_THROWS(to<char>("€"), std::logic_error);
    // Sanity check internal details
    test_utf_name();
    test_simple_encodings();
    test_win_codepages();

    test_latin1_conversions();
    test_utf_to_utf();
    test_utf_to_utf_allocator_support();

    std::cout << "Testing charset to/from UTF conversion functions\n";
    std::cout << "  char" << std::endl;
    test_utf_for<char>();
    std::cout << "  wchar_t" << std::endl;
    test_utf_for<wchar_t>();
#ifndef BOOST_LOCALE_NO_CXX20_STRING8
    std::cout << "  char8_t" << std::endl;
    test_utf_for<char8_t>();
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR16_T
    std::cout << "  char16_t" << std::endl;
    test_utf_for<char16_t>();
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR32_T
    std::cout << "  char32_t" << std::endl;
    test_utf_for<char32_t>();
#endif

    test_all_combinations();
    test_between();
}

// Internal tests, keep those out of the above scope

bool isLittleEndian()
{
#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__)
    return __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__;
#elif defined(__LITTLE_ENDIAN__)
    return true;
#elif defined(__BIG_ENDIAN__)
    return false;
#endif
    const int endianMark = 1;
    return reinterpret_cast<const char*>(&endianMark)[0] == 1;
}

#include "../src/boost/locale/util/encoding.hpp"
#include "../src/boost/locale/util/win_codepages.hpp"

void test_utf_name()
{
    TEST_EQ(boost::locale::util::utf_name<char>(), std::string("UTF-8"));
#ifdef __cpp_char8_t
    TEST_EQ(boost::locale::util::utf_name<char8_t>(), std::string("UTF-8"));
#endif
    TEST_EQ(boost::locale::util::utf_name<char16_t>(), std::string(isLittleEndian() ? "UTF-16LE" : "UTF-16BE"));
    TEST_EQ(boost::locale::util::utf_name<char32_t>(), std::string(isLittleEndian() ? "UTF-32LE" : "UTF-32BE"));
}

void test_simple_encodings()
{
    using namespace boost::locale::util;
    const auto encodings = get_simple_encodings();
    for(auto it = encodings.begin(), end = encodings.end(); it != end; ++it) {
        TEST_EQ(normalize_encoding(*it), *it); // Must be normalized
        const auto it2 = std::find(it + 1, end, *it);
        TEST(it2 == end);
        if(it2 != end)
            std::cerr << "Duplicate entry: " << *it << '\n'; // LCOV_EXCL_LINE
    }
    const auto it = std::is_sorted_until(encodings.begin(), encodings.end());
    TEST(it == encodings.end());
    if(it != encodings.end())
        std::cerr << "First wrongly sorted element: " << *it << '\n'; // LCOV_EXCL_LINE
}

void test_win_codepages()
{
    using namespace boost::locale::util;

    for(const windows_encoding *it = all_windows_encodings, *end = std::end(all_windows_encodings); it != end; ++it) {
        TEST_EQ(normalize_encoding(it->name), it->name); // Must be normalized
        auto is_same_win_codepage = [&it](const windows_encoding& rhs) -> bool {
            return it->codepage == rhs.codepage && std::strcmp(it->name, rhs.name) == 0;
        };
        const auto* it2 = std::find_if(it + 1, end, is_same_win_codepage);
        TEST(it2 == end);
        if(it2 != end)
            std::cerr << "Duplicate entry: " << it->name << ':' << it->codepage << '\n'; // LCOV_EXCL_LINE
    }
    const auto cmp = [](const windows_encoding& rhs, const windows_encoding& lhs) -> bool { return rhs < lhs.name; };
    const auto* it = std::is_sorted_until(all_windows_encodings, std::end(all_windows_encodings), cmp);
    TEST(it == std::end(all_windows_encodings));
    if(it != std::end(all_windows_encodings))
        std::cerr << "First wrongly sorted element: " << it->name << '\n'; // LCOV_EXCL_LINE
}

// boostinspect:noascii
