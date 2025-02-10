//
// Copyright (c) 2022-2023 Alexander Grund
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/locale/hold_ptr.hpp>
#include <boost/locale/util.hpp>
#include <boost/locale/util/locale_data.hpp>
#include "boostLocale/test/test_helpers.hpp"
#include "boostLocale/test/tools.hpp"
#include "boostLocale/test/unit_test.hpp"
#include <cstdlib>
#include <stdexcept>

namespace {
struct Dummy {
    int i_;
    Dummy(int i) : i_(i) { ++ctr; }
    ~Dummy() { --ctr; }
    Dummy(const Dummy&) = delete;
    Dummy(Dummy&&) = delete;

    int foo() { return i_; }
    int foo() const { return -i_; }

    static int ctr;
};
int Dummy::ctr = 0;
} // namespace

void test_hold_ptr()
{
    {
        boost::locale::hold_ptr<Dummy> empty;
        TEST(!empty);
        auto* raw = new Dummy(42);
        boost::locale::hold_ptr<Dummy> ptr(raw);
        const boost::locale::hold_ptr<Dummy>& const_ptr = ptr;
        TEST_REQUIRE(ptr);
        TEST(ptr.get() == raw);
        TEST(const_ptr.get() == raw);
        // const propagation
        TEST_EQ((*ptr).foo(), raw->i_);
        TEST_EQ((*const_ptr).foo(), -raw->i_);
        TEST_EQ(ptr->foo(), raw->i_);
        TEST_EQ(const_ptr->foo(), -raw->i_);
        TEST_EQ(ptr.get()->foo(), raw->i_);
        TEST_EQ(const_ptr.get()->foo(), -raw->i_);
        // move construct
        boost::locale::hold_ptr<Dummy> ptr2 = std::move(ptr);
        TEST(!ptr);
        TEST_REQUIRE(ptr2);
        TEST(ptr2.get() == raw);
        // move assign
        ptr = std::move(ptr2);
        TEST(ptr);
        TEST_REQUIRE(!ptr2);
        TEST(ptr.get() == raw);
        // Swap
        boost::locale::hold_ptr<Dummy> ptr3(new Dummy(1337));
        ptr.swap(ptr3);
        TEST_EQ(ptr->foo(), 1337);
        TEST_EQ(ptr3->foo(), 42);
    }
    TEST_EQ(Dummy::ctr, 0);
    auto* raw = new Dummy(42);
    {
        boost::locale::hold_ptr<Dummy> ptr(new Dummy(1));
        TEST_EQ(Dummy::ctr, 2);
        ptr.reset(raw);
        TEST_EQ(Dummy::ctr, 1);
        TEST_EQ(ptr->foo(), 42);
        TEST(ptr.release() == raw);
        TEST_EQ(Dummy::ctr, 1);
    }
    TEST_EQ(Dummy::ctr, 1);
    delete raw;
}

void test_get_system_locale()
{
    // Clear all -> Default to C
    {
        using boost::locale::test::unsetenv;
        unsetenv("LC_CTYPE");
        unsetenv("LC_ALL");
        unsetenv("LANG");
    }

    using boost::locale::util::get_system_locale;
#if !BOOST_LOCALE_USE_WIN32_API
    TEST_EQ(get_system_locale(false), "C");
#else
    // On Windows the user default name is used, so we can only test the encoding
    TEST(get_system_locale(true).find(".UTF-8") != std::string::npos);
    {
        const std::string loc = get_system_locale(false);
        const std::string enc = loc.substr(loc.find_last_of('.'));
        // encoding should be a windows codepage, but in the error case can be UTF-8
        if(enc.find(".windows-") != 0u)
            TEST_EQ(enc, ".UTF-8"); // LCOV_EXCL_LINE
    }
#endif
    // LC_ALL, LC_CTYPE and LANG variables used in this order
    using boost::locale::test::setenv;
    setenv("LANG", "mylang.foo");
    TEST_EQ(get_system_locale(true), "mylang.foo");
    setenv("LC_CTYPE", "this.lang");
    TEST_EQ(get_system_locale(true), "this.lang");
    setenv("LC_ALL", "barlang.bar");
    TEST_EQ(get_system_locale(true), "barlang.bar");
}

void test_locale_data()
{
    boost::locale::util::locale_data data;
    // Default is C.US-ASCII
    TEST_EQ(data.language(), "C");
    TEST_EQ(data.country(), "");
    TEST_EQ(data.encoding(), "US-ASCII");
    TEST(!data.is_utf8());
    TEST_EQ(data.variant(), "");

    TEST(data.parse("en_US.UTF-8"));
    TEST_EQ(data.language(), "en");
    TEST_EQ(data.country(), "US");
    TEST_EQ(data.encoding(), "UTF-8");
    TEST(data.is_utf8());
    TEST_EQ(data.variant(), "");

    TEST(data.parse("C"));
    TEST_EQ(data.language(), "C");
    TEST_EQ(data.country(), "");
    TEST_EQ(data.encoding(), "US-ASCII");
    TEST(!data.is_utf8());
    TEST_EQ(data.variant(), "");

    TEST(data.parse("ku_TR.UTF-8@sorani"));
    TEST_EQ(data.language(), "ku");
    TEST_EQ(data.country(), "TR");
    TEST_EQ(data.encoding(), "UTF-8");
    TEST(data.is_utf8());
    TEST_EQ(data.variant(), "sorani");

    TEST(data.parse("POSIX"));
    TEST_EQ(data.language(), "C");
    TEST_EQ(data.country(), "");
    TEST_EQ(data.encoding(), "US-ASCII");
    TEST(!data.is_utf8());
    TEST_EQ(data.variant(), "");

    TEST(data.parse("da_DK.ISO8859-15@euro"));
    TEST_EQ(data.language(), "da");
    TEST_EQ(data.country(), "DK");
    TEST_EQ(data.encoding(), "ISO8859-15");
    TEST(!data.is_utf8());
    TEST_EQ(data.variant(), "euro");

    TEST(data.parse("de_DE.ISO8859-1"));
    TEST_EQ(data.language(), "de");
    TEST_EQ(data.country(), "DE");
    TEST_EQ(data.encoding(), "ISO8859-1");
    TEST(!data.is_utf8());
    TEST_EQ(data.variant(), "");

    TEST(data.parse("ja_JP.eucJP"));
    TEST_EQ(data.language(), "ja");
    TEST_EQ(data.country(), "JP");
    TEST_EQ(data.encoding(), "EUCJP");
    TEST(!data.is_utf8());
    TEST_EQ(data.variant(), "");

    TEST(data.parse("ko_KR.EUC@dict"));
    TEST_EQ(data.language(), "ko");
    TEST_EQ(data.country(), "KR");
    TEST_EQ(data.encoding(), "EUC");
    TEST(!data.is_utf8());
    TEST_EQ(data.variant(), "dict");

    TEST(data.parse("th_TH.TIS620"));
    TEST_EQ(data.language(), "th");
    TEST_EQ(data.country(), "TH");
    TEST_EQ(data.encoding(), "TIS620");
    TEST(!data.is_utf8());
    TEST_EQ(data.variant(), "");

    TEST(data.parse("zh_TW.UTF-8@radical"));
    TEST_EQ(data.language(), "zh");
    TEST_EQ(data.country(), "TW");
    TEST_EQ(data.encoding(), "UTF-8");
    TEST(data.is_utf8());
    TEST_EQ(data.variant(), "radical");

    // Country can be a 3-digit value
    TEST(data.parse("en_001.UTF-8"));
    TEST_EQ(data.language(), "en");
    TEST_EQ(data.country(), "001");
    TEST_EQ(data.encoding(), "UTF-8");
    TEST(data.is_utf8());
    TEST_EQ(data.variant(), "");

    // to_string yields the input (if format is correct already)
    for(const std::string name : {"C",
                                  "en_US.UTF-8",
                                  "ku_TR.UTF-8@sorani",
                                  "da_DK.ISO8859-15@euro",
                                  "de_DE.ISO8859-1",
                                  "en_US",
                                  "ko_KR.EUC@dict",
                                  "th_TH.TIS620",
                                  "zh_TW.UTF-8@radical",
                                  "en_001",
                                  "en_150.UTF-8"})
    {
        TEST(data.parse(name));
        TEST_EQ(data.to_string(), name);
    }
    // US-ASCII encoding is ignored
    TEST(data.parse("da_TR.US-ASCII"));
    TEST_EQ(data.to_string(), "da_TR");
    TEST(data.parse("da_TR.US-ASCII@dic"));
    TEST_EQ(data.to_string(), "da_TR@dic");

    // Unify casing:
    // - language: lowercase
    // - region: uppercase
    // - encoding: uppercase
    // - variant: lowercase
    TEST(data.parse("EN_us.utf-8@EUro"));
    TEST_EQ(data.language(), "en");
    TEST_EQ(data.country(), "US");
    TEST_EQ(data.encoding(), "UTF-8");
    TEST(data.is_utf8());
    TEST_EQ(data.variant(), "euro");
    TEST_EQ(data.to_string(), "en_US.UTF-8@euro");
    TEST(data.parse("lAnGUagE_cOunTRy.eNCo-d123inG@Va-r1_Ant"));
    TEST_EQ(data.to_string(), "language_COUNTRY.ENCO-D123ING@va-r1_ant");

    // Dash is allowed in addition to underscore
    TEST(data.parse("de-DE.UTF-8"));
    TEST_EQ(data.to_string(), "de_DE.UTF-8");

    // C/POSIX is allowed to have an encoding
    TEST(data.parse("C.UTF-8"));
    TEST_EQ(data.to_string(), "C.UTF-8");
    TEST(data.parse("POSIX.UTF-8"));
    TEST_EQ(data.to_string(), "C.UTF-8");

    // Special case: en_US_POSIX is an alias for "C"
    TEST(data.parse("en_US_POSIX"));
    TEST_EQ(data.to_string(), "C");
    TEST(data.parse("En_Us_POsix.UTF-8"));
    TEST_EQ(data.to_string(), "C.UTF-8");

    // Missing values are defaulted
    TEST(data.parse("en"));
    TEST_EQ(data.to_string(), "en");
    TEST_EQ(data.encoding(), "US-ASCII");
    TEST(!data.is_utf8());
    TEST(data.parse("en.UTF-8"));
    TEST_EQ(data.to_string(), "en.UTF-8");
    TEST_EQ(data.encoding(), "UTF-8");
    TEST(data.is_utf8());
    TEST(data.parse("en@dict"));
    TEST_EQ(data.to_string(), "en@dict");
    TEST_EQ(data.encoding(), "US-ASCII");
    TEST_EQ(data.variant(), "dict");
    TEST(data.parse("en_US@dict"));
    TEST_EQ(data.to_string(), "en_US@dict");
    TEST_EQ(data.encoding(), "US-ASCII");
    TEST_EQ(data.variant(), "dict");

    // Error cases, default values used starting from error

    // Invalid language (separator at start or not an ASCII letter)
    for(const std::string invalidName :
        {"_en_US.UTF-8", "-en_US.UTF-8", ".en_US.UTF-8", "@en_US.UTF-8", "e1_US.UTF-8", "eö_US.UTF-8"})
    {
        TEST(!data.parse(invalidName));
        TEST_EQ(data.to_string(), "C");
    }
    // Invalid country
    TEST(!data.parse("en_UÖ.UTF-8"));
    TEST_EQ(data.to_string(), "en");
    TEST(!data.parse("en_1234.UTF-8")); // To many digits
    TEST_EQ(data.to_string(), "en");
    TEST(!data.parse("en_US1.UTF-8")); // digits in text
    TEST_EQ(data.to_string(), "en");
    TEST(!data.parse("en_1US.UTF-8")); // digits in text
    TEST_EQ(data.to_string(), "en");

    // Empty parts:
    // Language
    TEST(!data.parse("_US.UTF-8@variant"));
    TEST_EQ(data.to_string(), "C");
    // Country
    TEST(!data.parse("en_.UTF-8@variant"));
    TEST_EQ(data.to_string(), "en");
    // Encoding
    TEST(!data.parse("en_US.@variant"));
    TEST_EQ(data.to_string(), "en_US");
    // Variant
    TEST(!data.parse("en_US.UTF-8@"));
    TEST_EQ(data.to_string(), "en_US.UTF-8");

    // C/POSIX with any other field except the encoding
    for(const std::string invalidName : {"C_US", "C@variant", "POSIX_US", "POSIX@variant"}) {
        TEST(!data.parse(invalidName));
        TEST_EQ(data.to_string(), "C");
    }

    // Construct from string
    TEST_EQ(boost::locale::util::locale_data("en_US.UTF-8").to_string(), "en_US.UTF-8");
    TEST_THROWS(boost::locale::util::locale_data invalid("en_UÖ.UTF-8"), std::invalid_argument);
}

#include "../src/boost/locale/util/numeric.hpp"
#include <limits>
#include <locale>
#include <sstream>

void test_try_to_int()
{
    using boost::locale::util::try_to_int;

    int v = 1337;
    TEST(try_to_int("0", v));
    TEST_EQ(v, 0);

    TEST(try_to_int("42", v));
    TEST_EQ(v, 42);

    TEST(try_to_int("-1337", v));
    TEST_EQ(v, -1337);

    std::ostringstream ss;
    ss.imbue(std::locale::classic());
    empty_stream(ss) << std::numeric_limits<int>::min();
    TEST(try_to_int(ss.str(), v));
    TEST_EQ(v, std::numeric_limits<int>::min());
    empty_stream(ss) << std::numeric_limits<int>::max();
    TEST(try_to_int(ss.str(), v));
    TEST_EQ(v, std::numeric_limits<int>::max());

    TEST(!try_to_int("", v));
    TEST(!try_to_int("a", v));
    TEST(!try_to_int("1.", v));
    TEST(!try_to_int("1a", v));
    TEST(!try_to_int("a1", v));
    static_assert(sizeof(long long) > sizeof(int), "Value below under/overflows!");
    empty_stream(ss) << static_cast<long long>(std::numeric_limits<int>::min()) - 1;
    TEST(!try_to_int(ss.str(), v));
    empty_stream(ss) << static_cast<long long>(std::numeric_limits<int>::max()) + 1;
    TEST(!try_to_int(ss.str(), v));
}

void test_main(int /*argc*/, char** /*argv*/)
{
    test_hold_ptr();
    test_get_system_locale();
    test_locale_data();
    test_try_to_int();
}
