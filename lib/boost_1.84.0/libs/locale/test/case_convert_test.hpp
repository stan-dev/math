//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
// Copyright (c) 2022-2023 Alexander Grund
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/locale/conversion.hpp>
#include "boostLocale/test/tools.hpp"
#include <locale>
#include <string>

namespace boost { namespace locale { namespace case_convert_test {

    template<typename CharType>
    void test_one_impl(const std::locale& l,
                       const std::basic_string<CharType>& src,
                       const std::basic_string<CharType>& tgt_lower,
                       const std::basic_string<CharType>& tgt_upper)
    {
        TEST_EQ(boost::locale::normalize(src, boost::locale::norm_default, l), src); // Already normalized
        TEST_EQ(boost::locale::to_upper(src, l), tgt_upper);
        TEST_EQ(boost::locale::to_lower(src, l), tgt_lower);
        TEST_EQ(boost::locale::fold_case(src, l), tgt_lower);
    }

    template<typename CharType>
    void
    test_one(const std::locale& l, const std::string& src, const std::string& tgt_lower, const std::string& tgt_upper)
    {
        test_one_impl(l,
                      to_correct_string<CharType>(src, l),
                      to_correct_string<CharType>(tgt_lower, l),
                      to_correct_string<CharType>(tgt_upper, l));
    }

    template<typename CharType>
    void test_no_op_title_case(const std::locale& l, const std::string& src)
    {
        const auto src_c = to_correct_string<CharType>(src, l);
        const auto src_l = boost::locale::to_lower(src, l);
        const auto src_u = boost::locale::to_upper(src, l);
        TEST_EQ(boost::locale::to_title(src_c, l), src_c);
        TEST_EQ(boost::locale::to_title(src_l, l), src_l);
        TEST_EQ(boost::locale::to_title(src_u, l), src_u);
    }

}}} // namespace boost::locale::case_convert_test
