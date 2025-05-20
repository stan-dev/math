//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
// Copyright (c) 2021-2023 Alexander Grund
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/locale/encoding.hpp>
#include <boost/locale/generator.hpp>
#include <boost/locale/gnu_gettext.hpp>
#include <boost/locale/localization_backend.hpp>
#include <boost/locale/message.hpp>
#include "boostLocale/test/tools.hpp"
#include "boostLocale/test/unit_test.hpp"
#include <fstream>
#include <iostream>
#include <limits>
#include <type_traits>
#include <vector>

namespace bl = boost::locale;

void test_messages_info()
{
    using string_vec = std::vector<std::string>;
    {
        bl::gnu_gettext::messages_info info;
        info.locale_category = "LC";
        TEST_EQ(info.get_catalog_paths(), string_vec{});
        info.paths.push_back(".");
        TEST_EQ(info.get_catalog_paths(), string_vec{"./C/LC"});
        info.language = "en";
        TEST_EQ(info.get_catalog_paths(), string_vec{"./en/LC"});
        info.country = "US";
        TEST_EQ(info.get_catalog_paths(), (string_vec{"./en_US/LC", "./en/LC"}));
        info.country.clear();
        info.variant = "euro";
        TEST_EQ(info.get_catalog_paths(), (string_vec{"./en@euro/LC", "./en/LC"}));
        info.country = "US";
        TEST_EQ(info.get_catalog_paths(), (string_vec{"./en_US@euro/LC", "./en@euro/LC", "./en_US/LC", "./en/LC"}));

        info.paths = string_vec{"/1", "/2"};
        TEST_EQ(info.get_catalog_paths(),
                (string_vec{"/1/en_US@euro/LC",
                            "/2/en_US@euro/LC",
                            "/1/en@euro/LC",
                            "/2/en@euro/LC",
                            "/1/en_US/LC",
                            "/2/en_US/LC",
                            "/1/en/LC",
                            "/2/en/LC"}));
    }
}

std::string backend;
bool file_loader_is_actually_called = false;

struct file_loader {
    std::vector<char> operator()(const std::string& name, const std::string& /*encoding*/) const
    {
        std::vector<char> buffer;
        std::ifstream f(name.c_str(), std::ifstream::binary);
        if(!f)
            return buffer;
        f.seekg(0, std::ifstream::end);
        const auto len = f.tellg();
        f.seekg(0);
        if(len > 0) {
            buffer.resize(static_cast<size_t>(len));
            f.read(buffer.data(), buffer.size());
        }
        file_loader_is_actually_called = true;
        return buffer;
    }
};

std::string same_s(std::string s)
{
    return s;
}

std::wstring same_w(std::wstring s)
{
    return s;
}

#ifndef BOOST_LOCALE_NO_CXX20_STRING8
std::basic_string<char8_t> same_u8(std::basic_string<char8_t> s)
{
    return s;
}
#endif

#ifdef BOOST_LOCALE_ENABLE_CHAR16_T
std::u16string same_u16(std::u16string s)
{
    return s;
}
#endif

#ifdef BOOST_LOCALE_ENABLE_CHAR32_T
std::u32string same_u32(std::u32string s)
{
    return s;
}
#endif

namespace impl {

template<class T, std::size_t = sizeof(T)>
std::true_type is_complete_impl(T*);
std::false_type is_complete_impl(...);

template<class T>
using has_ctype = decltype(is_complete_impl(std::declval<std::ctype<T>*>()));

template<typename Char, typename... Ts>
typename std::enable_if<has_ctype<Char>::value, bool>::type stream_translate(std::basic_ostream<Char>& ss, Ts&&... args)
{
    ss << bl::translate(args...);
    return true;
}

// Required for char types not fully supported by the standard library
// e.g.: error: implicit instantiation of undefined template 'std::ctype<char8_t>'
template<typename Char, typename... Ts>
typename std::enable_if<!has_ctype<Char>::value, bool>::type stream_translate(std::basic_ostream<Char>&, Ts&&...)
{
    return false; // LCOV_EXCL_LINE
}

template<typename Char>
void test_cntranslate(const std::string& sContext,
                      const std::string& sSingular,
                      const std::string& sPlural,
                      long long n,
                      const std::string& sExpected,
                      const std::locale& l,
                      const std::string& domain)
{
    typedef std::basic_string<Char> string_type;
    const string_type expected = to_correct_string<Char>(sExpected, l);

    const string_type c = to<Char>(sContext);
    const string_type s = to<Char>(sSingular);
    const string_type p = to<Char>(sPlural);

    if(domain == "default") {
        TEST_EQ(bl::translate(c, s, p, n).str(l), expected);
        TEST_EQ(bl::translate(c.c_str(), s.c_str(), p.c_str(), n).str(l), expected);
        std::locale tmp_locale;
        std::locale::global(l);
        string_type tmp = bl::translate(c, s, p, n);
        TEST_EQ(tmp, expected);
        tmp = bl::translate(c, s, p, n).str();
        TEST_EQ(tmp, expected);
        std::locale::global(tmp_locale);

        std::basic_ostringstream<Char> ss;
        ss.imbue(l);
        if(stream_translate(ss, c, s, p, n))
            TEST_EQ(ss.str(), expected);

        // Copyable & movable
        const string_type s2 = ascii_to<Char>("missing Singular");
        const string_type p2 = ascii_to<Char>("missing Plural");
        const string_type& expected2 = (n == 1) ? s2 : p2;
        auto translation1 = bl::translate(c, s, p, n);
        auto translation2 = bl::translate(c, s2, p2, n);
        TEST_EQ(translation1.str(l), expected);
        TEST_EQ(translation2.str(l), expected2);
        // Copy
        translation1 = translation2;
        TEST_EQ(translation1.str(l), expected2);
        {
            bl::basic_message<Char> t3(translation2);
            TEST_EQ(t3.str(l), expected2);
        }
        // Move
        translation1 = bl::translate(c, s, p, n);
        translation2 = std::move(translation1);
        TEST_EQ(translation2.str(l), expected);
        {
            bl::basic_message<Char> t3(std::move(translation2));
            TEST_EQ(t3.str(l), expected);
        }
        // Swap
        translation1 = bl::translate(c, s, p, n);
        translation2 = bl::translate(c, s2, p2, n);
        TEST_EQ(translation1.str(l), expected);
        TEST_EQ(translation2.str(l), expected2);
        using std::swap;
        swap(translation1, translation2);
        TEST_EQ(translation1.str(l), expected2);
        TEST_EQ(translation2.str(l), expected);
        translation1 = bl::translate(c, s, p, n);
        translation2 = bl::translate(c, s2, p2, 1); // n==1!
        swap(translation1, translation2);
        TEST_EQ(translation1.str(l), s2);
        TEST_EQ(translation2.str(l), expected);
    }
    TEST_EQ(bl::translate(c, s, p, n).str(l, domain), expected);
    std::locale tmp_locale;
    std::locale::global(l);
    TEST_EQ(bl::translate(c, s, p, n).str(domain), expected);
    std::locale::global(tmp_locale);
    {
        std::basic_ostringstream<Char> ss;
        ss.imbue(l);
        if(stream_translate(ss << bl::as::domain(domain), c, s, p, n))
            TEST_EQ(ss.str(), expected);
    }
    {
        std::basic_ostringstream<Char> ss;
        ss.imbue(l);
        if(stream_translate(ss << bl::as::domain(domain), c.c_str(), s.c_str(), p.c_str(), n))
            TEST_EQ(ss.str(), expected);
    }
    // Missing facet -> No translation
    {
        const string_type nonAscii = ((string_type() + Char('\x7F')) + Char('\x82')) + Char('\xF0');
        const string_type p2 = nonAscii + p + nonAscii;
        // For char the non-ASCII chars are removed -> original p
        const string_type expected2 = (n == 1) ? s : (std::is_same<Char, char>::value ? p : p2);
        TEST_EQ(bl::translate(c, s, p2, n).str(std::locale::classic()), expected2);
    }
}

template<typename Char>
void test_ntranslate(const std::string& sSingular,
                     const std::string& sPlural,
                     long long n,
                     const std::string& sExpected,
                     const std::locale& l,
                     const std::string& domain)
{
    typedef std::basic_string<Char> string_type;
    const string_type expected = to_correct_string<Char>(sExpected, l);
    const string_type s = to<Char>(sSingular);
    const string_type p = to<Char>(sPlural);
    if(domain == "default") {
        TEST_EQ(bl::translate(s, p, n).str(l), expected);
        TEST_EQ(bl::translate(s.c_str(), p.c_str(), n).str(l), expected);
        std::locale tmp_locale;
        std::locale::global(l);
        string_type tmp = bl::translate(s, p, n);
        TEST_EQ(tmp, expected);
        tmp = bl::translate(s, p, n).str();
        TEST_EQ(tmp, expected);
        std::locale::global(tmp_locale);

        std::basic_ostringstream<Char> ss;
        ss.imbue(l);
        if(stream_translate(ss, s, p, n))
            TEST_EQ(ss.str(), expected);
    }
    TEST_EQ(bl::translate(s, p, n).str(l, domain), expected);
    std::locale tmp_locale;
    std::locale::global(l);
    TEST_EQ(bl::translate(s, p, n).str(domain), expected);
    std::locale::global(tmp_locale);
    {
        std::basic_ostringstream<Char> ss;
        ss.imbue(l);
        if(stream_translate(ss << bl::as::domain(domain), s, p, n))
            TEST_EQ(ss.str(), expected);
    }
    {
        std::basic_ostringstream<Char> ss;
        ss.imbue(l);
        if(stream_translate(ss << bl::as::domain(domain), s.c_str(), p.c_str(), n))
            TEST_EQ(ss.str(), expected);
    }
}

template<typename Char>
void test_ctranslate(const std::string& sContext,
                     const std::string& sOriginal,
                     const std::string& sExpected,
                     const std::locale& l,
                     const std::string& domain)
{
    typedef std::basic_string<Char> string_type;
    const string_type expected = to_correct_string<Char>(sExpected, l);
    const string_type original = to<Char>(sOriginal);
    const string_type c = to<Char>(sContext);
    if(domain == "default") {
        TEST_EQ(bl::translate(c, original).str(l), expected);
        TEST_EQ(bl::translate(c.c_str(), original.c_str()).str(l), expected);
        std::locale tmp_locale;
        std::locale::global(l);
        string_type tmp = bl::translate(c, original);
        TEST_EQ(tmp, expected);
        tmp = bl::translate(c, original).str();
        TEST_EQ(tmp, expected);
        std::locale::global(tmp_locale);

        std::basic_ostringstream<Char> ss;
        ss.imbue(l);
        if(stream_translate(ss, c, original))
            TEST_EQ(ss.str(), expected);
    }
    TEST_EQ(bl::translate(c, original).str(l, domain), expected);
    std::locale tmp_locale;
    std::locale::global(l);
    TEST_EQ(bl::translate(c, original).str(domain), expected);
    std::locale::global(tmp_locale);
    {
        std::basic_ostringstream<Char> ss;
        ss.imbue(l);
        if(stream_translate(ss << bl::as::domain(domain), c, original))
            TEST_EQ(ss.str(), expected);
    }
    {
        std::basic_ostringstream<Char> ss;
        ss.imbue(l);
        if(stream_translate(ss << bl::as::domain(domain), c.c_str(), original.c_str()))
            TEST_EQ(ss.str(), expected);
    }
}

template<typename Char>
void test_translate(const std::string& sOriginal,
                    const std::string& sExpected,
                    const std::locale& l,
                    const std::string& domain)
{
    typedef std::basic_string<Char> string_type;
    const string_type expected = to_correct_string<Char>(sExpected, l);
    const string_type original = to<Char>(sOriginal);
    if(domain == "default") {
        TEST_EQ(bl::translate(original).str(l), expected);
        TEST_EQ(bl::translate(original.c_str()).str(l), expected);
        std::locale tmp_locale;
        std::locale::global(l);
        string_type tmp = bl::translate(original);
        TEST_EQ(tmp, expected);
        tmp = bl::translate(original).str();
        TEST_EQ(tmp, expected);
        std::locale::global(tmp_locale);

        std::basic_ostringstream<Char> ss;
        ss.imbue(l);
        if(stream_translate(ss, original))
            TEST_EQ(ss.str(), expected);
    }
    TEST_EQ(bl::translate(original).str(l, domain), expected);
    std::locale tmp_locale;
    std::locale::global(l);
    TEST_EQ(bl::translate(original).str(domain), expected);
    std::locale::global(tmp_locale);
    {
        std::basic_ostringstream<Char> ss;
        ss.imbue(l);
        if(stream_translate(ss << bl::as::domain(domain), original))
            TEST_EQ(ss.str(), expected);
    }
    {
        std::basic_ostringstream<Char> ss;
        ss.imbue(l);
        if(stream_translate(ss << bl::as::domain(domain), original.c_str()))
            TEST_EQ(ss.str(), expected);
    }
}
} // namespace impl

void test_cntranslate(const std::string& c,
                      const std::string& s,
                      const std::string& p,
                      long long n,
                      const std::string& expected,
                      const std::locale& l,
                      const std::string& domain)
{
    std::cout << "  char" << std::endl;
    impl::test_cntranslate<char>(c, s, p, n, expected, l, domain);
    std::cout << "  wchar_t" << std::endl;
    impl::test_cntranslate<wchar_t>(c, s, p, n, expected, l, domain);
#ifndef BOOST_LOCALE_NO_CXX20_STRING8
    std::cout << "  char8_t" << std::endl;
    impl::test_cntranslate<char8_t>(c, s, p, n, expected, l, domain);
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR16_T
    std::cout << "  char16_t" << std::endl;
    impl::test_cntranslate<char16_t>(c, s, p, n, expected, l, domain);
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR32_T
    std::cout << "  char32_t" << std::endl;
    impl::test_cntranslate<char32_t>(c, s, p, n, expected, l, domain);
#endif
}

void test_ntranslate(const std::string& s,
                     const std::string& p,
                     long long n,
                     const std::string& expected,
                     const std::locale& l,
                     const std::string& domain)
{
    std::cout << "  char" << std::endl;
    impl::test_ntranslate<char>(s, p, n, expected, l, domain);
    std::cout << "  wchar_t" << std::endl;
    impl::test_ntranslate<wchar_t>(s, p, n, expected, l, domain);
#ifndef BOOST_LOCALE_NO_CXX20_STRING8
    std::cout << "  char8_t" << std::endl;
    impl::test_ntranslate<char8_t>(s, p, n, expected, l, domain);
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR16_T
    std::cout << "  char16_t" << std::endl;
    impl::test_ntranslate<char16_t>(s, p, n, expected, l, domain);
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR32_T
    std::cout << "  char32_t" << std::endl;
    impl::test_ntranslate<char32_t>(s, p, n, expected, l, domain);
#endif
}

void test_ctranslate(const std::string& c,
                     const std::string& original,
                     const std::string& expected,
                     const std::locale& l,
                     const std::string& domain)
{
    std::cout << "  char" << std::endl;
    impl::test_ctranslate<char>(c, original, expected, l, domain);
    std::cout << "  wchar_t" << std::endl;
    impl::test_ctranslate<wchar_t>(c, original, expected, l, domain);
#ifndef BOOST_LOCALE_NO_CXX20_STRING8
    std::cout << "  char8_t" << std::endl;
    impl::test_ctranslate<char8_t>(c, original, expected, l, domain);
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR16_T
    std::cout << "  char16_t" << std::endl;
    impl::test_ctranslate<char16_t>(c, original, expected, l, domain);
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR32_T
    std::cout << "  char32_t" << std::endl;
    impl::test_ctranslate<char32_t>(c, original, expected, l, domain);
#endif
}

void test_translate(const std::string& original,
                    const std::string& expected,
                    const std::locale& l,
                    const std::string& domain)
{
    std::cout << "  char" << std::endl;
    impl::test_translate<char>(original, expected, l, domain);
    std::cout << "  wchar_t" << std::endl;
    impl::test_translate<wchar_t>(original, expected, l, domain);
#ifndef BOOST_LOCALE_NO_CXX20_STRING8
    std::cout << "  char8_t" << std::endl;
    impl::test_translate<char8_t>(original, expected, l, domain);
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR16_T
    std::cout << "  char16_t" << std::endl;
    impl::test_translate<char16_t>(original, expected, l, domain);
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR32_T
    std::cout << "  char32_t" << std::endl;
    impl::test_translate<char32_t>(original, expected, l, domain);
#endif
}

bool iso_8859_8_supported = true;

void test_main(int argc, char** argv)
{
    test_messages_info();

    const std::string message_path = (argc == 2) ? argv[1] : ".";

    for(const std::string& backend_name : boost::locale::localization_backend_manager::global().get_all_backends()) {
        std::cout << "Testing for backend --------- " << backend_name << std::endl;
        boost::locale::localization_backend_manager tmp_backend = boost::locale::localization_backend_manager::global();
        tmp_backend.select(backend_name);
        boost::locale::localization_backend_manager::global(tmp_backend);

        backend = backend_name;

        boost::locale::generator g;
        g.add_messages_domain("simple");
        g.add_messages_domain("full");
        // Fallback using only ASCII keys, choose a specific encoding != UTF-8, here: Latin1
        g.add_messages_domain("fall/ISO-8859-1");
        g.add_messages_path(message_path);
        g.set_default_messages_domain("default");

        for(const std::string locale_name : {"he_IL.UTF-8", "he_IL.ISO8859-8"}) {
            std::locale l;

            if(locale_name.find(".ISO") != std::string::npos) {
                try {
                    l = g(locale_name);
                } catch(const boost::locale::conv::invalid_charset_error&) {                      // LCOV_EXCL_LINE
                    std::cout << "Looks like ISO-8859-8 is not supported! skipping" << std::endl; // LCOV_EXCL_LINE
                    iso_8859_8_supported = false;                                                 // LCOV_EXCL_LINE
                    continue;                                                                     // LCOV_EXCL_LINE
                }
            } else
                l = g(locale_name);

            std::cout << "  Testing " << locale_name << std::endl;
            std::cout << "    single forms" << std::endl;

            test_translate("hello", "שלום", l, "default");
            test_translate("hello", "היי", l, "simple");
            test_translate("hello", "hello", l, "undefined");
            test_translate("untranslated", "untranslated", l, "default");
            // Check removal of old "context" information
            test_translate("#untranslated", "#untranslated", l, "default");
            test_translate("##untranslated", "##untranslated", l, "default");
            test_ctranslate("context", "hello", "שלום בהקשר אחר", l, "default");
            test_translate("#hello", "#שלום", l, "default");

            std::cout << "    plural forms" << std::endl;

            {
                test_ntranslate("x day", "x days", 0, "x ימים", l, "default");
                test_ntranslate("x day", "x days", 1, "יום x", l, "default");
                test_ntranslate("x day", "x days", 2, "יומיים", l, "default");
                test_ntranslate("x day", "x days", 3, "x ימים", l, "default");
                test_ntranslate("x day", "x days", 20, "x יום", l, "default");
                test_ntranslate("x day", "x days", 0, "x days", l, "undefined");
                test_ntranslate("x day", "x days", 1, "x day", l, "undefined");
                test_ntranslate("x day", "x days", 2, "x days", l, "undefined");
                test_ntranslate("x day", "x days", 20, "x days", l, "undefined");
                // Ensure no truncation occurs
                test_ntranslate("x day", "x days", std::numeric_limits<long long>::min(), "x days", l, "undefined");
                test_ntranslate("x day", "x days", std::numeric_limits<long long>::max(), "x days", l, "undefined");
                for(unsigned bit = 1; bit < std::numeric_limits<long long>::digits; ++bit) {
                    // Set each individual bit possible and add 1.
                    // If the value is truncated the 1 will remain leading to singular form
                    const auto num = static_cast<long long>(static_cast<unsigned long long>(1) << bit);
                    test_ntranslate("x day", "x days", num + 1, "x days", l, "undefined");
                }
            }
            std::cout << "    plural forms with context" << std::endl;
            {
                std::string inp = "context";
                std::string out = "בהקשר ";

                test_cntranslate(inp, "x day", "x days", 0, out + "x ימים", l, "default");
                test_cntranslate(inp, "x day", "x days", 1, out + "יום x", l, "default");
                test_cntranslate(inp, "x day", "x days", 2, out + "יומיים", l, "default");
                test_cntranslate(inp, "x day", "x days", 3, out + "x ימים", l, "default");
                test_cntranslate(inp, "x day", "x days", 20, out + "x יום", l, "default");

                test_cntranslate(inp, "x day", "x days", 0, "x days", l, "undefined");
                test_cntranslate(inp, "x day", "x days", 1, "x day", l, "undefined");
                test_cntranslate(inp, "x day", "x days", 2, "x days", l, "undefined");
                test_cntranslate(inp, "x day", "x days", 20, "x days", l, "undefined");
            }
        }
        std::cout << "  Testing fallbacks" << std::endl;
        {
            const std::locale l = g("he_IL.UTF-8");
            test_translate("test", "he_IL", l, "full");
            test_translate("test", "he", l, "fall");
            for(int n = -1; n < 5; ++n) {
                // No plural forms -> Use english logic
                // Singular is translated, plural is not
                test_ntranslate("test", "tests", n, (n == 1) ? "he" : "tests", l, "fall");
            }
        }

        std::cout << "  Testing automatic conversions " << std::endl;
        std::locale::global(g("he_IL.UTF-8"));

        TEST_EQ(same_s(bl::translate("hello")), "שלום");
        TEST_EQ(same_w(bl::translate(to<wchar_t>("hello"))), to<wchar_t>("שלום"));

#ifndef BOOST_LOCALE_NO_CXX20_STRING8
        TEST_EQ(same_u8(bl::translate(to<char8_t>("hello"))), to<char8_t>("שלום"));
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR16_T
        if(backend == "icu" || backend == "std")
            TEST_EQ(same_u16(bl::translate(to<char16_t>("hello"))), to<char16_t>("שלום"));
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR32_T
        if(backend == "icu" || backend == "std")
            TEST_EQ(same_u32(bl::translate(to<char32_t>("hello"))), to<char32_t>("שלום"));
#endif
    }

    std::cout << "Testing custom file system support" << std::endl;
    {
        boost::locale::gnu_gettext::messages_info info;
        info.language = "he";
        info.country = "IL";
        info.encoding = "UTF-8";
        info.paths.push_back(message_path);

        info.domains.push_back(bl::gnu_gettext::messages_info::domain("default"));
        info.callback = file_loader();

        file_loader_is_actually_called = false;
        std::locale l(std::locale::classic(), boost::locale::gnu_gettext::create_messages_facet<char>(info));
        TEST(file_loader_is_actually_called);
        TEST_EQ(bl::translate("hello").str(l), "שלום");
    }
    if(iso_8859_8_supported) {
        std::cout << "Testing non-US-ASCII keys" << std::endl;
        std::cout << "  UTF-8 keys" << std::endl;
        {
            boost::locale::generator g;
            g.add_messages_domain("default");
            g.add_messages_path(message_path);

            std::locale l = g("he_IL.UTF-8");

            // narrow
            TEST_EQ(bl::gettext("בדיקה", l), "test");
            TEST_EQ(bl::gettext("לא קיים", l), "לא קיים");

            // wide
            std::wstring wtest = bl::conv::utf_to_utf<wchar_t>("בדיקה");
            std::wstring wmiss = bl::conv::utf_to_utf<wchar_t>("לא קיים");
            TEST_EQ(bl::gettext(wtest.c_str(), l), L"test");
            TEST_EQ(bl::gettext(wmiss.c_str(), l), wmiss);

            l = g("he_IL.ISO-8859-8");

            // conversion with substitution
            TEST_EQ(bl::gettext("test-あにま-בדיקה", l), bl::conv::from_utf("test--בדיקה", "ISO-8859-8"));
        }

        std::cout << "  ANSI keys" << std::endl;

        {
            boost::locale::generator g;
            g.add_messages_domain("default/ISO-8859-8");
            g.add_messages_path(message_path);

            std::locale l = g("he_IL.UTF-8");

            // narrow non-UTF-8 keys
            // match
            TEST_EQ(bl::gettext(bl::conv::from_utf("בדיקה", "ISO-8859-8").c_str(), l), "test");
            // conversion
            TEST_EQ(bl::gettext(bl::conv::from_utf("לא קיים", "ISO-8859-8").c_str(), l), "לא קיים");
        }
    }
    // Test compiles
    {
        bl::gettext("");
        bl::gettext(L"");
        bl::dgettext("", "");
        bl::dgettext("", L"");

        bl::pgettext("", "");
        bl::pgettext(L"", L"");
        bl::dpgettext("", "", "");
        bl::dpgettext("", L"", L"");

        bl::ngettext("", "", 1);
        bl::ngettext(L"", L"", 1);
        bl::dngettext("", "", "", 1);
        bl::dngettext("", L"", L"", 1);

        bl::npgettext("", "", "", 1);
        bl::npgettext(L"", L"", L"", 1);
        bl::dnpgettext("", "", "", "", 1);
        bl::dnpgettext("", L"", L"", L"", 1);
    }
}

// boostinspect:noascii
