// Copyright 2024 Jeff Trull.
//
// Distributed under the Boost Software License, Version 1.0.
//
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt

// There are some tricky corner cases that occur when Wave handles large inputs.
// Tokens at, near, or crossing the boundary at which a "fill" occurs can
// interact in unexpected ways. This test is intended to increase code
// coverage for such situations.

#include <boost/wave/cpplexer/re2clex/cpp_re.hpp>   // for BOOST_WAVE_BSIZE

#include <boost/wave.hpp>
#include <boost/wave/token_ids.hpp>
#include <boost/wave/cpplexer/cpp_lex_token.hpp>
#include <boost/wave/cpplexer/cpp_lex_iterator.hpp>

#include <string>

using token_t = boost::wave::cpplexer::lex_token<>;
using lex_iter_t = boost::wave::cpplexer::lex_iterator<token_t>;
using ctx_t = boost::wave::context<
    std::string::iterator, lex_iter_t,
    boost::wave::iteration_context_policies::load_file_to_string,
    boost::wave::context_policies::default_preprocessing_hooks>;


int main() {
    using namespace boost::wave;

    // construct input with a trigraph in a tricky spot
    constexpr std::size_t space_count = BOOST_WAVE_BSIZE - 18;
    auto inp_txt = std::string(space_count, ' ') + ";";
    // add in pound trigraph one character at a time
    // otherwise, the string literal may get translated by the compiler
    inp_txt.push_back('?');
    inp_txt.push_back('?');
    inp_txt.push_back('=');
    inp_txt += "                  auto foo = bar;\n";

    ctx_t ctx( inp_txt.begin(), inp_txt.end(), "longfile.cpp" );

    auto it = ctx.begin();

    if (token_id(*it) != T_SPACE)
        return 1;

    if (it->get_value().size() != space_count)
        return 2;

    ++it;

    if (token_id(*it) != T_SEMICOLON)
        return 3;

    ++it;

    if (token_id(*it) != T_POUND_TRIGRAPH)
        return 4;

    ++it;

    if (token_id(*it) != T_SPACE)
        return 5;

    ++it;

    if (token_id(*it) != T_AUTO)
        return 6;

    return 0;
}
