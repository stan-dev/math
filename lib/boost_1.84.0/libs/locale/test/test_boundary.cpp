//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#define BOOST_LOCALE_ERROR_LIMIT 100000

#include <boost/locale/boundary.hpp>
#include <boost/locale/generator.hpp>
#include <boost/locale/localization_backend.hpp>
#include "boostLocale/test/tools.hpp"
#include "boostLocale/test/unit_test.hpp"
#include <boost/assert.hpp>
#include <iostream>
#include <list>
#include <vector>
#ifdef BOOST_LOCALE_WITH_ICU
#    include <unicode/uversion.h>
#endif

namespace lb = boost::locale::boundary;
template<typename Char>
using chunks_t = std::vector<std::basic_string<Char>>;
using masks_t = std::vector<unsigned>;
using positions_t = std::vector<size_t>;

template<typename Iterator, typename Char>
void run_segment_iterator_test(const lb::segment_index<Iterator>& map,
                               const Iterator begin,
                               const Iterator end,
                               const chunks_t<Char>& chunks,
                               const masks_t& masks,
                               const positions_t& pos)
{
    {
        unsigned i = 0;
        typename lb::segment_index<Iterator>::iterator p;
        for(p = map.begin(); p != map.end(); ++p, i++) {
            TEST_REQUIRE(i < masks.size());
            TEST_EQ(p->str(), chunks[i]);
            TEST_EQ(p->rule(), masks[i]);
        }
        TEST_EQ(i, chunks.size());

        for(;;) {
            if(p == map.begin()) {
                TEST_EQ(i, 0u);
                break;
            } else {
                --p, --i;
                TEST_EQ(p->str(), chunks[i]);
                TEST_EQ(p->rule(), masks[i]);
            }
        }
        for(i = 0, p = map.end(); i < chunks.size(); i++) {
            --p;
            size_t index = chunks.size() - i - 1;
            TEST_EQ(p->str(), chunks[index]);
            TEST_EQ(p->rule(), masks[index]);
        }
        TEST(p == map.begin());
    }

    {
        size_t chunk_ptr = 0, i = 0;
        for(Iterator optr = begin; optr != end; optr++, i++) {
            const auto p = map.find(optr);
            if(chunk_ptr < pos.size() && i >= pos[chunk_ptr])
                chunk_ptr++;
            if(chunk_ptr >= pos.size())
                TEST(p == map.end());
            else {
                TEST_EQ(p->str(), chunks[chunk_ptr]);
                TEST_EQ(p->rule(), unsigned(masks[chunk_ptr]));
            }
        }
    }
}

template<typename Iterator>
void run_break_iterator_test(const lb::boundary_point_index<Iterator>& map,
                             const Iterator begin,
                             const Iterator end,
                             const std::vector<Iterator>& iters,
                             const masks_t& masks)
{
    unsigned i = 0;
    typename lb::boundary_point_index<Iterator>::iterator p;
    for(p = map.begin(); p != map.end(); ++p, i++) {
        TEST_REQUIRE(i < masks.size());
        TEST(p->iterator() == iters[i]);
        TEST_EQ(p->rule(), masks[i]);
    }

    TEST_EQ(i, iters.size());

    do {
        --p;
        --i;
        TEST(p->iterator() == iters.at(i));
    } while(p != map.begin());
    TEST_EQ(i, 0u);

    unsigned iters_ptr = 0;
    for(Iterator optr = begin; optr != end; optr++) {
        p = map.find(optr);
        TEST(p->iterator() == iters[iters_ptr]);
        if(iters.at(iters_ptr) == optr)
            iters_ptr++;
    }
}

template<typename Iterator>
void verify_index(const lb::boundary_point_index<Iterator>& map,
                  const std::vector<Iterator>& iters,
                  const masks_t& masks)
{
    BOOST_ASSERT(iters.size() == masks.size());
    TEST_REQUIRE(static_cast<size_t>(std::distance(map.begin(), map.end())) == masks.size());
    size_t i = 0;
    for(const auto& b_point : map) {
        TEST(b_point.iterator() == iters[i]);
        TEST_EQ(b_point.rule(), masks[i]);
        ++i;
    }
}

template<typename Iterator, typename Char>
void verify_index(const lb::segment_index<Iterator>& map, const chunks_t<Char>& chunks, const masks_t& masks)
{
    BOOST_ASSERT(chunks.size() == masks.size());
    TEST_REQUIRE(static_cast<size_t>(std::distance(map.begin(), map.end())) == masks.size());
    size_t i = 0;
    for(const auto& seg : map) {
        TEST_EQ(seg.str(), chunks[i]);
        TEST_EQ(seg.rule(), masks[i]);
        ++i;
    }
}

template<typename Char, typename Iterator>
void test_word_container(Iterator begin,
                         Iterator end,
                         const std::vector<size_t>& ipos,
                         const std::vector<unsigned>& imasks,
                         const std::vector<std::basic_string<Char>>& ichunks,
                         std::locale l,
                         lb::boundary_type bt = lb::word)
{
    using segments_t = lb::segment_index<Iterator>;
    using boundaries_t = lb::boundary_point_index<Iterator>;
    for(int sm = (bt == lb::word ? 31 : 3); sm >= 0; sm--) {
        unsigned mask = ((sm & 1) != 0) * 0xF + ((sm & 2) != 0) * 0xF0 + ((sm & 4) != 0) * 0xF00
                        + ((sm & 8) != 0) * 0xF000 + ((sm & 16) != 0) * 0xF0000;

        masks_t masks;
        std::vector<size_t> pos;
        std::vector<unsigned> boundary_masks;
        std::basic_string<Char> empty_chunk;

        chunks_t<Char> chunks;
        chunks_t<Char> full_chunks;
        std::vector<Iterator> iters;
        iters.push_back(begin);
        boundary_masks.push_back(0);

        for(unsigned i = 0; i < imasks.size(); i++) {
            if(imasks[i] & mask) {
                masks.push_back(imasks[i]);
                chunks.push_back(ichunks[i]);
                full_chunks.push_back(empty_chunk + ichunks[i]);
                empty_chunk.clear();
                pos.push_back(ipos[i]);
            } else
                empty_chunk += ichunks[i];

            if((imasks[i] & mask) || i == imasks.size() - 1) {
                Iterator ptr = begin;
                std::advance(ptr, ipos[i]);
                iters.push_back(ptr);
                boundary_masks.push_back(imasks[i]);
            }
        }
        {
            segments_t map(bt, begin, end, l);
            map.rule(mask);
            map.full_select(false);
            run_segment_iterator_test(map, begin, end, chunks, masks, pos);
            map.full_select(true);
            run_segment_iterator_test(map, begin, end, full_chunks, masks, pos);
        }
        {
            boundaries_t map(bt, begin, end, l);
            map.rule(mask);
            run_break_iterator_test(map, begin, end, iters, boundary_masks);
        }

        std::cout << "-- Copy from segment_index\n";
        {
            segments_t ti(bt, begin, end, l);
            ti.rule(mask);
            std::cout << "---- Construct boundary_point_index\n";
            {
                boundaries_t bi(ti);
                bi.rule(mask);
                verify_index(bi, iters, boundary_masks);
            }
            std::cout << "---- Assign boundary_point_index\n";
            {
                boundaries_t bi;
                bi.rule(mask);
                bi = ti;
                verify_index(bi, iters, boundary_masks);
            }
            std::cout << "---- Construct segment_index\n";
            {
                segments_t ti2(ti);
                verify_index(ti2, chunks, masks);
            }
            std::cout << "---- Assign segment_index\n";
            {
                segments_t ti2;
                ti2 = ti;
                verify_index(ti2, chunks, masks);
            }
        }
        std::cout << "-- Copy from boundary_point_index\n";
        {
            boundaries_t bi(bt, begin, end, l);
            bi.rule(mask);
            std::cout << "---- Construct boundary_point_index\n";
            {
                boundaries_t bi2(bi);
                verify_index(bi2, iters, boundary_masks);
            }
            std::cout << "---- Assign boundary_point_index\n";
            {
                boundaries_t bi2;
                bi2 = bi;
                verify_index(bi2, iters, boundary_masks);
            }
            std::cout << "---- Construct segment_index\n";
            {
                segments_t ti(bi);
                ti.rule(mask);
                verify_index(ti, chunks, masks);
            }
            std::cout << "---- Assign segment_index\n";
            {
                segments_t ti;
                ti.rule(mask);
                ti = bi;
                verify_index(ti, chunks, masks);
            }
        }
    }
}

template<typename Char>
void run_word(std::string* original,
              const int* none,
              const int* num,
              const int* word,
              const int* kana,
              const int* ideo,
              std::locale l,
              lb::boundary_type b = lb::word)
{
    std::vector<size_t> pos;
    std::vector<std::basic_string<Char>> chunks;
    std::vector<unsigned> masks;
    std::basic_string<Char> test_string;
    for(int i = 0; !original[i].empty(); i++) {
        chunks.push_back(to_correct_string<Char>(original[i], l));
        test_string += chunks.back();
        pos.push_back(test_string.size());
        masks.push_back((none && none[i] ? 0xFu : 0u) | (num && num[i] ? 0xF0u : 0u) | (word && word[i] ? 0xF00u : 0u)
                        | (kana && kana[i] ? 0xF000u : 0u) | (ideo && ideo[i] ? 0xF0000u : 0u));
    }

    std::list<Char> lst(test_string.begin(), test_string.end());
    test_word_container<Char>(lst.begin(), lst.end(), pos, masks, chunks, l, b);
    test_word_container<Char>(test_string.begin(), test_string.end(), pos, masks, chunks, l, b);
}

std::string character[] = {"שָ", "ל", "וֹ", "ם", "!", ""};
int nones[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

// clang-format off
std::string sentence1[]={"To be\n","or not\n","to be?\n"," That is the question. ","Or maybe not",""};
int         sentence1a[]={      0,          0,        1,                         1,             0, 0};
int         sentence1b[]={      1,          1,        0,                         0,             1, 0};

std::string line1[]={"To ","be\n","or ","not\n","to ","be",""};
int         line1a[]={ 1,   0,     1 ,  0,       1,   1 , 0 };
int         line1b[]={ 0,   1,     0 ,  1,       0,   0 , 0 };
// clang-format on

void test_boundaries(std::string* all, int* first, int* second, lb::boundary_type t)
{
    boost::locale::generator g;
    std::cout << " char UTF-8" << std::endl;
    run_word<char>(all, first, second, nullptr, nullptr, nullptr, g("he_IL.UTF-8"), t);
    std::cout << " char CP1255" << std::endl;
    run_word<char>(all, first, second, nullptr, nullptr, nullptr, g("he_IL.cp1255"), t);
    std::cout << " wchar_t" << std::endl;
    run_word<wchar_t>(all, first, second, nullptr, nullptr, nullptr, g("he_IL.UTF-8"), t);
#ifndef BOOST_LOCALE_NO_CXX20_STRING8
    std::cout << " char8_t" << std::endl;
    run_word<char8_t>(all, first, second, nullptr, nullptr, nullptr, g("he_IL.UTF-8"), t);
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR16_T
    std::cout << " char16_t" << std::endl;
    run_word<char16_t>(all, first, second, nullptr, nullptr, nullptr, g("he_IL.UTF-8"), t);
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR32_T
    std::cout << " char32_t" << std::endl;
    run_word<char32_t>(all, first, second, nullptr, nullptr, nullptr, g("he_IL.UTF-8"), t);
#endif
}

void word_boundary()
{
    boost::locale::generator g;
    int zero[25] = {0};
    std::string txt_empty[] = {""};

    // clang-format off
    std::string txt_simple[] = {" ","Hello",",","World","!"," ",""};
    int        none_simple[] = { 1,      0,  1,      0,  1,  1, 0};
    int        word_simple[] = { 0,      1,  0,      1,  0,  0, 0};

    std::string txt_all[] = {"10"," ","Hello"," ","Windows7"," ","He22o"," ","平仮名","アヒル",""};
    int        none_all[] = {  0,  1,      0,  1,         0,  1,      0,  1,      0,      0,  0};
#if U_ICU_VERSION_MAJOR_NUM >= 62
    // ICU 62+ returns only the number classification if there is a number at the boundary
    int         num_all[] = {  1,  0,      0,  0,         1,  0,      0,  0,      0,      0,  0};
    int        word_all[] = {  0,  0,      1,  0,         0,  0,      1,  0,      0,      0,  0};
#else
    // ICU < 62 combines the word and number classification if there is a number at the boundary
    int         num_all[] = {  1,  0,      0,  0,         1,  0,      0,  0,      0,      0,  0}; // LCOV_EXCL_LINE
    int        word_all[] = {  0,  0,      1,  0,         1,  0,      1,  0,      0,      0,  0}; // LCOV_EXCL_LINE
#endif
#if U_ICU_VERSION_MAJOR_NUM >= 50
    int        kana_all[] = {  0,  0,      0,  0,         0,  0,      0,  0,      0,      0,  0};
    int        ideo_all[] = {  0,  0,      0,  0,         0,  0,      0,  0,      1,      1,  1};
#else
    int        kana_all[] = {  0,  0,      0,  0,         0,  0,      0,  0,      0,      1,  1}; // LCOV_EXCL_LINE
    int        ideo_all[] = {  0,  0,      0,  0,         0,  0,      0,  0,      1,      0,  0}; // LCOV_EXCL_LINE
#endif
    // clang-format on

    std::cout << " char UTF-8" << std::endl;
    const std::locale utf8_en_locale = g("en_US.UTF-8");
    const std::locale utf8_jp_locale = g("ja_JP.UTF-8");
    run_word<char>(txt_empty, zero, zero, zero, zero, zero, utf8_en_locale);
    run_word<char>(txt_simple, none_simple, zero, word_simple, zero, zero, utf8_en_locale);
    run_word<char>(txt_all, none_all, num_all, word_all, kana_all, ideo_all, utf8_jp_locale);

    std::cout << " char Shift-JIS" << std::endl;
    const std::locale sjis_jp_locale = g("ja_JP.SJIS");
    run_word<char>(txt_empty, zero, zero, zero, zero, zero, sjis_jp_locale);
    run_word<char>(txt_simple, none_simple, zero, word_simple, zero, zero, sjis_jp_locale);
    run_word<char>(txt_all, none_all, num_all, word_all, kana_all, ideo_all, sjis_jp_locale);

    std::cout << " wchar_t" << std::endl;
    run_word<wchar_t>(txt_empty, zero, zero, zero, zero, zero, utf8_en_locale);
    run_word<wchar_t>(txt_simple, none_simple, zero, word_simple, zero, zero, utf8_en_locale);
    run_word<wchar_t>(txt_all, none_all, num_all, word_all, kana_all, ideo_all, utf8_jp_locale);

#ifndef BOOST_LOCALE_NO_CXX20_STRING8
    std::cout << " char8_t" << std::endl;
    run_word<char8_t>(txt_empty, zero, zero, zero, zero, zero, g("ja_JP.UTF-8"));
    run_word<char8_t>(txt_simple, none_simple, zero, word_simple, zero, zero, utf8_en_locale);
    run_word<char8_t>(txt_all, none_all, num_all, word_all, kana_all, ideo_all, utf8_jp_locale);
#endif

#ifdef BOOST_LOCALE_ENABLE_CHAR16_T
    std::cout << " char16_t" << std::endl;
    run_word<char16_t>(txt_empty, zero, zero, zero, zero, zero, g("ja_JP.UTF-8"));
    run_word<char16_t>(txt_simple, none_simple, zero, word_simple, zero, zero, utf8_en_locale);
    run_word<char16_t>(txt_all, none_all, num_all, word_all, kana_all, ideo_all, utf8_jp_locale);
#endif

#ifdef BOOST_LOCALE_ENABLE_CHAR32_T
    std::cout << " char32_t" << std::endl;
    run_word<char32_t>(txt_empty, zero, zero, zero, zero, zero, g("ja_JP.UTF-8"));
    run_word<char32_t>(txt_simple, none_simple, zero, word_simple, zero, zero, utf8_en_locale);
    run_word<char32_t>(txt_all, none_all, num_all, word_all, kana_all, ideo_all, utf8_jp_locale);
#endif
}
void test_op_one_side(const std::string& sl, const std::string& sr, int val)
{
    boost::locale::boundary::ssegment l(sl.begin(), sl.end(), 0), r(sr.begin(), sr.end(), 0);

#if BOOST_LOCALE_SPACESHIP_NULLPTR_WARNING
#    pragma clang diagnostic push
#    pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif

    // segment
    TEST_EQ((l == r), (val == 0));
    TEST_EQ((l != r), (val != 0));
    TEST_EQ((l <= r), (val <= 0));
    TEST_EQ((l < r), (val < 0));
    TEST_EQ((l >= r), (val >= 0));
    TEST_EQ((l > r), (val > 0));

    // C string
    TEST_EQ((l == sr.c_str()), (val == 0));
    TEST_EQ((l != sr.c_str()), (val != 0));
    TEST_EQ((l <= sr.c_str()), (val <= 0));
    TEST_EQ((l < sr.c_str()), (val < 0));
    TEST_EQ((l >= sr.c_str()), (val >= 0));
    TEST_EQ((l > sr.c_str()), (val > 0));

    TEST_EQ((sl.c_str() == r), (val == 0));
    TEST_EQ((sl.c_str() != r), (val != 0));
    TEST_EQ((sl.c_str() <= r), (val <= 0));
    TEST_EQ((sl.c_str() < r), (val < 0));
    TEST_EQ((sl.c_str() >= r), (val >= 0));
    TEST_EQ((sl.c_str() > r), (val > 0));

    // C++ string
    TEST_EQ((l == sr), (val == 0));
    TEST_EQ((l != sr), (val != 0));
    TEST_EQ((l <= sr), (val <= 0));
    TEST_EQ((l < sr), (val < 0));
    TEST_EQ((l >= sr), (val >= 0));
    TEST_EQ((l > sr), (val > 0));

    TEST_EQ((sl == r), (val == 0));
    TEST_EQ((sl != r), (val != 0));
    TEST_EQ((sl <= r), (val <= 0));
    TEST_EQ((sl < r), (val < 0));
    TEST_EQ((sl >= r), (val >= 0));
    TEST_EQ((sl > r), (val > 0));
    // self check
    TEST_EQ((sl == sr), (val == 0));
    TEST_EQ((sl != sr), (val != 0));
    TEST_EQ((sl <= sr), (val <= 0));
    TEST_EQ((sl < sr), (val < 0));
    TEST_EQ((sl >= sr), (val >= 0));
    TEST_EQ((sl > sr), (val > 0));

#if BOOST_LOCALE_SPACESHIP_NULLPTR_WARNING
#    pragma clang diagnostic pop
#endif
}

void test_op(const std::string& sl, const std::string& sr, int val)
{
    test_op_one_side(sl, sr, val);
    test_op_one_side(sr, sl, -val);
}
void segment_operator()
{
    test_op("", "a", -1);
    test_op("", "", 0);
    test_op("aa", "aaa", -1);
    test_op("aa", "ab", -1);
}

BOOST_LOCALE_DISABLE_UNREACHABLE_CODE_WARNING
void test_main(int /*argc*/, char** /*argv*/)
{
#ifndef BOOST_LOCALE_NO_STD_BACKEND
    {
        namespace bl = boost::locale;
        const bl::localization_backend_manager orig_backend = bl::localization_backend_manager::global();
        bl::localization_backend_manager tmp_backend = bl::localization_backend_manager::global();
        tmp_backend.select("std");
        bl::localization_backend_manager::global(tmp_backend);

        bl::generator g;
        const std::string text = "To be or not to be, that is the question.";
        // Std backend doesn't support segmentation, expect reasonable error
        TEST_THROWS(bl::boundary::ssegment_index map(bl::boundary::word, text.begin(), text.end(), g("en_US.UTF-8")),
                    std::runtime_error);
        bl::localization_backend_manager::global(orig_backend);
    }
#endif
#ifndef BOOST_LOCALE_WITH_ICU
    std::cout << "ICU is not build... Skipping\n";
    return;
#endif // !BOOST_LOCALE_WITH_ICU
    std::cout << "Testing segment operators" << std::endl;
    segment_operator();
    std::cout << "Testing word boundary" << std::endl;
    word_boundary();
    std::cout << "Testing character boundary" << std::endl;
    test_boundaries(character, nones, nullptr, lb::character);
    std::cout << "Testing sentence boundary" << std::endl;
    test_boundaries(sentence1, sentence1a, sentence1b, lb::sentence);
    std::cout << "Testing line boundary" << std::endl;
    test_boundaries(line1, line1a, line1b, lb::line);
}

// boostinspect:noascii
