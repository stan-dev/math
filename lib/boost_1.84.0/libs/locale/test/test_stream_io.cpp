//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
// Copyright (c) 2022-2023 Alexander Grund
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/locale/generator.hpp>
#include <boost/locale/localization_backend.hpp>
#include <fstream>

#include "boostLocale/test/tools.hpp"
#include "boostLocale/test/unit_test.hpp"

bool test_iso;
const bool test_iso_8859_8 =
#if defined(BOOST_LOCALE_WITH_ICU) || defined(BOOST_LOCALE_WITH_ICONV)
  true;
#else
  hasWinCodepage(28598);
#endif
bool test_utf;
bool test_sjis;

std::string he_il_8bit;
std::string en_us_8bit;
std::string ja_jp_shiftjis;

template<typename Char>
std::basic_string<Char> read_file(std::basic_istream<Char>& in)
{
    std::basic_string<Char> res;
    Char c;
    while(in.get(c))
        res += c;
    return res;
}

template<typename Char>
void test_ok(const std::string& content, const std::locale& l, std::basic_string<Char> cmp = std::basic_string<Char>())
{
    typedef std::basic_fstream<Char> stream_type;
    if(cmp.empty())
        cmp = to<Char>(content);

    {
        const std::string file_path = boost::locale::test::exe_name + "-test_read.txt";
        remove_file_on_exit _(file_path);
        {
            std::ofstream out_file(file_path, std::ios::binary);
            out_file << content;
        }
        stream_type in_file(file_path, stream_type::in);
        in_file.imbue(l);
        TEST_EQ(read_file<Char>(in_file), cmp);
    }

    {
        const std::string file_path = boost::locale::test::exe_name + "-test_write.txt";
        remove_file_on_exit _(file_path);
        {
            stream_type out_file(file_path, stream_type::out);
            out_file.imbue(l);
            out_file << cmp;
        }
        std::ifstream in_file(file_path);
        TEST_EQ(read_file<char>(in_file), content);
    }
}

template<typename Char>
void test_read_fail(const std::string& content, const std::locale& l, const int pos)
{
    const std::string file_path = boost::locale::test::exe_name + "-test.txt";
    remove_file_on_exit _(file_path);
    {
        std::ofstream f(file_path, std::ios::binary);
        f << content;
    }

    std::basic_fstream<Char> f(file_path, std::ios::in);
    f.imbue(l);
    // Read up until the position
    for(int i = 0; i < pos && f; i++) {
        Char c;
        f.get(c);
    }
    // Reading may fail before the position, e.g. when the implementation reads more than the current char and hence
    // detects the error early
    if(f) {
        // Usually it succeeds so far and must fail now
        Char c;
        TEST(!f.get(c));
    }
}

template<typename Char>
void test_write_fail(const std::string& content, const std::locale& l, int pos)
{
    const std::string file_path = boost::locale::test::exe_name + "-test.txt";
    remove_file_on_exit _(file_path);

    typedef std::basic_fstream<Char> stream_type;
    stream_type f(file_path, stream_type::out);
    f.imbue(l);
    std::basic_string<Char> out = to<Char>(content);
    for(int i = 0; i < pos; i++) {
        f << out.at(i) << std::flush;
        TEST(f);
    }
    f << out.at(pos);
    TEST(f.fail() || (f << std::flush).fail());
}

template<typename Char>
void test_for_char()
{
    boost::locale::generator g;
    if(test_utf) {
        std::cout << "    UTF-8" << std::endl;
        test_ok<Char>("grüße\nn i", g("en_US.UTF-8"));
        test_read_fail<Char>("abc\xFF\xFF", g("en_US.UTF-8"), 3);
        std::cout << "    Testing codepoints above 0xFFFF" << std::endl;
        std::cout << "      Single U+2008A" << std::endl;
        test_ok<Char>("\xf0\xa0\x82\x8a", g("en_US.UTF-8")); // U+2008A
        std::cout << "      Single U+2008A within text" << std::endl;
        test_ok<Char>("abc\"\xf0\xa0\x82\x8a\"", g("en_US.UTF-8")); // U+2008A
        std::string one = "\xf0\xa0\x82\x8a";
        std::string res;
        for(unsigned i = 0; i < 1000; i++)
            res += one;
        std::cout << "      U+2008A x 1000" << std::endl;
        test_ok<Char>(res.c_str(), g("en_US.UTF-8")); // U+2008A
    }

    if(test_iso) {
        if(test_iso_8859_8) {
            std::cout << "    ISO8859-8" << std::endl;
            test_ok<Char>("hello \xf9\xec\xe5\xed", g(he_il_8bit), to<Char>("hello שלום"));
        }
        std::cout << "    ISO8859-1" << std::endl;
        test_ok<Char>(to<char>("grüße\nn i"), g(en_us_8bit), to<Char>("grüße\nn i"));
        test_write_fail<Char>("grüßen שלום", g(en_us_8bit), 7);
    }

    if(test_sjis) {
        std::cout << "    Shift-JIS" << std::endl;
        test_ok<Char>("\x93\xfa\x96\x7b",
                      g(ja_jp_shiftjis),
                      boost::locale::conv::utf_to_utf<Char>("\xe6\x97\xa5\xe6\x9c\xac")); // Japan
    }
}
void test_wide_io()
{
    std::cout << "  wchar_t" << std::endl;
    test_for_char<wchar_t>();

    // std::codecvt<char8_t doesn't have proper library support (e.g. MSVC doesn't export the id)

#if defined BOOST_LOCALE_ENABLE_CHAR16_T
    std::cout << "  char16_t" << std::endl;
    test_for_char<char16_t>();
#endif
#if defined BOOST_LOCALE_ENABLE_CHAR32_T
    std::cout << "  char32_t" << std::endl;
    test_for_char<char32_t>();
#endif
}

void test_main(int /*argc*/, char** /*argv*/)
{
    for(const std::string& backendName : boost::locale::localization_backend_manager::global().get_all_backends()) {
        boost::locale::localization_backend_manager tmp_backend = boost::locale::localization_backend_manager::global();
        tmp_backend.select(backendName);
        boost::locale::localization_backend_manager::global(tmp_backend);

        en_us_8bit = "en_US.ISO8859-1";
        he_il_8bit = "he_IL.ISO8859-8";
        ja_jp_shiftjis = "ja_JP.SJIS";
        if(backendName == "std") {
            en_us_8bit = get_std_name(en_us_8bit);
            he_il_8bit = get_std_name(he_il_8bit);
            std::string real_ja_jp_shiftjis;
            ja_jp_shiftjis = get_std_name(ja_jp_shiftjis, &real_ja_jp_shiftjis);
            if(!ja_jp_shiftjis.empty() && !test_std_supports_SJIS_codecvt(real_ja_jp_shiftjis))
                ja_jp_shiftjis = ""; // LCOV_EXCL_LINE
        }

        std::cout << "Testing for backend " << backendName << std::endl;

        if(backendName == "std") {
            test_iso = !he_il_8bit.empty() && !en_us_8bit.empty();
            test_sjis = !ja_jp_shiftjis.empty();
            test_utf = !get_std_name("en_US.UTF-8").empty() && !get_std_name("he_IL.UTF-8").empty();
        } else if(backendName == "winapi") {
            test_iso = false;
            test_sjis = false;
            test_utf = true;
        } else if(backendName == "posix") {
#ifdef BOOST_LOCALE_NO_POSIX_BACKEND
            throw std::logic_error("Unexpected backend"); // LCOV_EXCL_LINE
#else
            test_iso = has_posix_locale(he_il_8bit) && has_posix_locale(en_us_8bit);
            test_utf = has_posix_locale("en_US.UTF-8");
#    ifdef BOOST_LOCALE_WITH_ICONV
            test_sjis = has_posix_locale(ja_jp_shiftjis);
#    else
            test_sjis = false;
#    endif
#endif
        } else {
            test_iso = true;
            test_sjis = true;
            test_utf = true;
        }

        test_wide_io();
    }
}

// boostinspect:noascii
