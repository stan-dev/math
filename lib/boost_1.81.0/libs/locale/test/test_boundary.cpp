//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#define BOOST_LOCALE_ERROR_LIMIT 100000

#include <boost/locale/boundary.hpp>
#include <boost/locale/generator.hpp>
#include "boostLocale/test/tools.hpp"
#include "boostLocale/test/unit_test.hpp"
#include <iostream>
#include <list>
#include <vector>
#ifdef BOOST_LOCALE_WITH_ICU
#    include <unicode/uversion.h>
#endif

namespace lb = boost::locale::boundary;

template<typename Char, typename Iterator>
void test_word_container(Iterator begin,
                         Iterator end,
                         const std::vector<size_t>& ipos,
                         const std::vector<unsigned>& imasks,
                         const std::vector<std::basic_string<Char>>& ichunks,
                         std::locale l,
                         lb::boundary_type bt = lb::word)
{
    for(int sm = (bt == lb::word ? 31 : 3); sm >= 0; sm--) {
        unsigned mask = ((sm & 1) != 0) * 0xF + ((sm & 2) != 0) * 0xF0 + ((sm & 4) != 0) * 0xF00
                        + ((sm & 8) != 0) * 0xF000 + ((sm & 16) != 0) * 0xF0000;

        std::vector<unsigned> masks;
        std::vector<size_t> pos;
        std::vector<unsigned> bmasks;
        std::basic_string<Char> empty_chunk;

        std::vector<std::basic_string<Char>> chunks;
        std::vector<std::basic_string<Char>> fchunks;
        std::vector<Iterator> iters;
        iters.push_back(begin);
        bmasks.push_back(0);

        for(unsigned i = 0; i < imasks.size(); i++) {
            if(imasks[i] & mask) {
                masks.push_back(imasks[i]);
                chunks.push_back(ichunks[i]);
                fchunks.push_back(empty_chunk + ichunks[i]);
                empty_chunk.clear();
                pos.push_back(ipos[i]);
            } else {
                empty_chunk += ichunks[i];
            }

            if((imasks[i] & mask) || i == imasks.size() - 1) {
                Iterator ptr = begin;
                std::advance(ptr, ipos[i]);
                iters.push_back(ptr);
                bmasks.push_back(imasks[i]);
            }
        }

        // segment iterator tests
        {
            lb::segment_index<Iterator> map(bt, begin, end, l);
            typedef typename lb::segment_index<Iterator>::iterator iter_type;

            map.rule(mask);

            {
                unsigned i = 0;
                iter_type p;
                map.full_select(false);
                for(p = map.begin(); p != map.end(); ++p, i++) {
                    TEST_REQUIRE(i < chunks.size());
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
                    TEST_EQ(p->rule(), unsigned(masks[index]));
                }
                TEST(p == map.begin());
            }

            {
                unsigned i = 0;
                iter_type p;
                map.full_select(true);
                for(p = map.begin(); p != map.end(); ++p, i++) {
                    TEST_EQ(p->str(), fchunks[i]);
                    TEST_EQ(p->rule(), masks[i]);
                }

                TEST(chunks.size() == i);

                for(;;) {
                    if(p == map.begin()) {
                        TEST(i == 0);
                        break;
                    } else {
                        --p;
                        TEST_EQ(p->str(), fchunks[--i]);
                        TEST_EQ(p->rule(), masks[i]);
                    }
                }

                for(i = 0, p = map.end(); i < chunks.size(); i++) {
                    --p;
                    size_t index = chunks.size() - i - 1u;
                    TEST_EQ(p->str(), fchunks[index]);
                    TEST_EQ(p->rule(), unsigned(masks[index]));
                }
                TEST(p == map.begin());
            }

            {
                iter_type p;
                unsigned chunk_ptr = 0;
                unsigned i = 0;
                map.full_select(false);
                for(Iterator optr = begin; optr != end; optr++, i++) {
                    p = map.find(optr);
                    if(chunk_ptr < pos.size() && i >= unsigned(pos[chunk_ptr])) {
                        chunk_ptr++;
                    }
                    if(chunk_ptr >= pos.size()) {
                        TEST(p == map.end());
                    } else {
                        TEST_EQ(p->str(), chunks[chunk_ptr]);
                        TEST_EQ(p->rule(), unsigned(masks[chunk_ptr]));
                    }
                }
            }
            {
                iter_type p;
                unsigned chunk_ptr = 0;
                unsigned i = 0;
                map.full_select(true);
                for(Iterator optr = begin; optr != end; optr++, i++) {
                    p = map.find(optr);
                    if(chunk_ptr < pos.size() && i >= unsigned(pos[chunk_ptr])) {
                        chunk_ptr++;
                    }
                    if(chunk_ptr >= pos.size()) {
                        TEST(p == map.end());
                    } else {
                        TEST_EQ(p->str(), fchunks[chunk_ptr]);
                        TEST_EQ(p->rule(), unsigned(masks[chunk_ptr]));
                    }
                }
            }

        } // segment iterator tests

        { // break iterator tests
            lb::boundary_point_index<Iterator> map(bt, begin, end, l);
            typedef typename lb::boundary_point_index<Iterator>::iterator iter_type;

            map.rule(mask);

            unsigned i = 0;
            iter_type p;
            for(p = map.begin(); p != map.end(); ++p, i++) {
                TEST(p->iterator() == iters[i]);
                TEST_EQ(p->rule(), bmasks[i]);
            }

            TEST(iters.size() == i);

            do {
                --p;
                --i;
                TEST(p->iterator() == iters.at(i));
            } while(p != map.begin());
            TEST(i == 0);

            unsigned iters_ptr = 0;
            for(Iterator optr = begin; optr != end; optr++) {
                p = map.find(optr);
                TEST(p->iterator() == iters[iters_ptr]);
                if(iters.at(iters_ptr) == optr)
                    iters_ptr++;
            }

        } // break iterator tests

        { // copy test
            typedef lb::segment_index<Iterator> ti_type;
            typedef lb::boundary_point_index<Iterator> bi_type;
            { // segment to bound
                ti_type ti(bt, begin, end, l);
                ti.rule(mask);
                {
                    bi_type bi(ti);
                    bi.rule(mask);
                    unsigned i = 0;
                    typename bi_type::iterator p;
                    for(p = bi.begin(); p != bi.end(); ++p, i++) {
                        TEST(p->iterator() == iters[i]);
                        TEST_EQ(p->rule(), bmasks[i]);
                    }
                }
                {
                    bi_type bi;
                    bi.rule(mask);
                    bi = ti;
                    unsigned i = 0;
                    typename bi_type::iterator p;
                    for(p = bi.begin(); p != bi.end(); ++p, i++) {
                        TEST(p->iterator() == iters[i]);
                        TEST_EQ(p->rule(), bmasks[i]);
                    }
                }
                // boundary_point to bound
                bi_type bi_2(bt, begin, end, l);
                bi_2.rule(mask);
                {
                    bi_type bi(bi_2);
                    unsigned i = 0;
                    typename bi_type::iterator p;
                    for(p = bi.begin(); p != bi.end(); ++p, i++) {
                        TEST(p->iterator() == iters[i]);
                        TEST_EQ(p->rule(), bmasks[i]);
                    }
                }
                {
                    bi_type bi;
                    bi = bi_2;
                    unsigned i = 0;
                    typename bi_type::iterator p;
                    for(p = bi.begin(); p != bi.end(); ++p, i++) {
                        TEST(p->iterator() == iters[i]);
                        TEST_EQ(p->rule(), bmasks[i]);
                    }
                }
            }
            { // boundary_point to segment
                bi_type bi(bt, begin, end, l);
                {
                    ti_type ti(bi);
                    ti.rule(mask);
                    unsigned i = 0;
                    typename ti_type::iterator p;
                    for(p = ti.begin(); p != ti.end(); ++p, i++) {
                        TEST(p->str() == chunks[i]);
                        TEST_EQ(p->rule(), masks[i]);
                    }
                }
                {
                    ti_type ti;
                    ti.rule(mask);
                    ti = (bi);
                    unsigned i = 0;
                    typename ti_type::iterator p;
                    for(p = ti.begin(); p != ti.end(); ++p, i++) {
                        TEST_EQ(p->str(), chunks[i]);
                        TEST_EQ(p->rule(), masks[i]);
                    }
                }
                ti_type ti_2(bt, begin, end, l);
                ti_2.rule(mask);
                {
                    ti_type ti(ti_2);
                    unsigned i = 0;
                    typename ti_type::iterator p;
                    for(p = ti.begin(); p != ti.end(); ++p, i++) {
                        TEST_EQ(p->str(), chunks[i]);
                        TEST_EQ(p->rule(), masks[i]);
                    }
                }
                {
                    ti_type ti;
                    ti = (ti_2);
                    unsigned i = 0;
                    typename ti_type::iterator p;
                    for(p = ti.begin(); p != ti.end(); ++p, i++) {
                        TEST_EQ(p->str(), chunks[i]);
                        TEST_EQ(p->rule(), masks[i]);
                    }
                }
            }
        }
    } // for mask
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
    run_word<char>(all, first, second, 0, 0, 0, g("he_IL.UTF-8"), t);
    std::cout << " char CP1255" << std::endl;
    run_word<char>(all, first, second, 0, 0, 0, g("he_IL.cp1255"), t);
    std::cout << " wchar_t" << std::endl;
    run_word<wchar_t>(all, first, second, 0, 0, 0, g("he_IL.UTF-8"), t);
#ifdef BOOST_LOCALE_ENABLE_CHAR16_T
    std::cout << " char16_t" << std::endl;
    run_word<char16_t>(all, first, second, 0, 0, 0, g("he_IL.UTF-8"), t);
#endif
#ifdef BOOST_LOCALE_ENABLE_CHAR32_T
    std::cout << " char32_t" << std::endl;
    run_word<char32_t>(all, first, second, 0, 0, 0, g("he_IL.UTF-8"), t);
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
    int         num_all[] = {  1,  0,      0,  0,         1,  0,      0,  0,      0,      0,  0};
    int        word_all[] = {  0,  0,      1,  0,         1,  0,      1,  0,      0,      0,  0};
#endif
#if U_ICU_VERSION_MAJOR_NUM >= 50
    int        kana_all[] = {  0,  0,      0,  0,         0,  0,      0,  0,      0,      0,  0};
    int        ideo_all[] = {  0,  0,      0,  0,         0,  0,      0,  0,      1,      1,  1};
#else
    int        kana_all[] = {  0,  0,      0,  0,         0,  0,      0,  0,      0,      1,  1};
    int        ideo_all[] = {  0,  0,      0,  0,         0,  0,      0,  0,      1,      0,  0};
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

    // segment
    TEST((l == r) == (val == 0));
    TEST((l != r) == (val != 0));
    TEST((l <= r) == (val <= 0));
    TEST((l < r) == (val < 0));
    TEST((l >= r) == (val >= 0));
    TEST((l > r) == (val > 0));

    // C string
    TEST((l == sr.c_str()) == (val == 0));
    TEST((l != sr.c_str()) == (val != 0));
    TEST((l <= sr.c_str()) == (val <= 0));
    TEST((l < sr.c_str()) == (val < 0));
    TEST((l >= sr.c_str()) == (val >= 0));
    TEST((l > sr.c_str()) == (val > 0));

    TEST((sl.c_str() == r) == (val == 0));
    TEST((sl.c_str() != r) == (val != 0));
    TEST((sl.c_str() <= r) == (val <= 0));
    TEST((sl.c_str() < r) == (val < 0));
    TEST((sl.c_str() >= r) == (val >= 0));
    TEST((sl.c_str() > r) == (val > 0));

    // C++ string
    TEST((l == sr) == (val == 0));
    TEST((l != sr) == (val != 0));
    TEST((l <= sr) == (val <= 0));
    TEST((l < sr) == (val < 0));
    TEST((l >= sr) == (val >= 0));
    TEST((l > sr) == (val > 0));

    TEST((sl == r) == (val == 0));
    TEST((sl != r) == (val != 0));
    TEST((sl <= r) == (val <= 0));
    TEST((sl < r) == (val < 0));
    TEST((sl >= r) == (val >= 0));
    TEST((sl > r) == (val > 0));
    // self check
    TEST((sl == sr) == (val == 0));
    TEST((sl != sr) == (val != 0));
    TEST((sl <= sr) == (val <= 0));
    TEST((sl < sr) == (val < 0));
    TEST((sl >= sr) == (val >= 0));
    TEST((sl > sr) == (val > 0));
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
#ifndef BOOST_LOCALE_WITH_ICU
    std::cout << "ICU is not build... Skipping\n";
    return;
#endif // !BOOST_LOCALE_WITH_ICU
    std::cout << "Testing segment operators" << std::endl;
    segment_operator();
    std::cout << "Testing word boundary" << std::endl;
    word_boundary();
    std::cout << "Testing character boundary" << std::endl;
    test_boundaries(character, nones, 0, lb::character);
    std::cout << "Testing sentence boundary" << std::endl;
    test_boundaries(sentence1, sentence1a, sentence1b, lb::sentence);
    std::cout << "Testing line boundary" << std::endl;
    test_boundaries(line1, line1a, line1b, lb::line);
}

// boostinspect:noascii
