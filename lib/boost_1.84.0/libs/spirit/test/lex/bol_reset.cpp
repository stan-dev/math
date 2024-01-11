//  Copyright (c) 2001-2010 Hartmut Kaiser
//  Copyright (c) 2010 Mathias Gaunard
// 
//  Distributed under the Boost Software License, Version 1.0. (See accompanying 
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// This test makes sure that the BOL state (begin of line) is properly reset
// if a token matched at the beginning of a line is discarded using 
// lex::pass_fail.

#include <boost/spirit/include/lex_lexertl.hpp>

#include <boost/core/lightweight_test.hpp>
#include <boost/phoenix/operator/self.hpp>
#include <boost/phoenix/statement/sequence.hpp>

#include <sstream>

namespace spirit = boost::spirit;
namespace lex = spirit::lex;

typedef char const* content_iterator;
typedef lex::lexertl::token<content_iterator> token_type;

struct lexer
  : lex::lexer<lex::lexertl::actor_lexer<token_type> >
{
    lexer() : word("^[a-zA-Z0-9]+$", 1)
    {
        self =  word [ 
                    lex::_state = "O" 
                ]
            |   lex::token_def<>("!.*$") [(
                    lex::_state = "O"
                  , lex::_pass = lex::pass_flags::pass_ignore 
                )]
            |   lex::token_def<>('\n', 2) [ 
                    lex::_state = "O" 
                ] 
            ;
        
        self("O") = 
                lex::token_def<>(".") [(
                    lex::_state = "INITIAL"
                  , lex::_pass = lex::pass_flags::pass_fail 
                )]
            ;
    }
    
    lex::token_def<> word;
};

typedef lexer::iterator_type token_iterator;

int main()
{
    std::string const s = "!foo\nbar\n!baz";
    
    content_iterator begin = s.data();
    content_iterator end = s.data() + s.size();
    
    lexer l;
    token_iterator begin2 = l.begin(begin, end);
    token_iterator end2 = l.end();
    
    std::size_t test_data[] = { 2, 1, 2 };
    std::size_t const test_data_size = sizeof(test_data)/sizeof(test_data[0]);

    token_iterator it = begin2;
    std::size_t i = 0;
    for (/**/; it != end2 && i < test_data_size; ++it, ++i)
    {
        BOOST_TEST(it->id() == test_data[i]);
    }
    BOOST_TEST(it == end2);
    BOOST_TEST(i == test_data_size);

    return boost::report_errors();
}
