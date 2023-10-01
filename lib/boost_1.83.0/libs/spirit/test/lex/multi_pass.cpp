//  Copyright (c) 2001-2010 Hartmut Kaiser
//  Copyright (c) 2010 Mathias Gaunard
// 
//  Distributed under the Boost Software License, Version 1.0. (See accompanying 
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// problem in multi_pass showing up in conjunction with certain lexer usage patterns

#include <boost/spirit/include/support_multi_pass.hpp>
#include <boost/spirit/include/lex_lexertl.hpp>

#include <boost/core/lightweight_test.hpp>
#include <boost/phoenix/operator/self.hpp>
#include <boost/phoenix/statement/sequence.hpp>

#include <sstream>

namespace spirit = boost::spirit;
namespace lex = spirit::lex;

typedef char const* content_iterator;

struct identifier
{
    identifier(content_iterator, content_iterator)
    {
    }
};

struct string_literal
{
    string_literal(content_iterator, content_iterator)
    {
    }
};

typedef lex::lexertl::token<
    content_iterator, boost::mpl::vector<identifier, string_literal>
> token_type;

struct lexer
  : lex::lexer<lex::lexertl::actor_lexer<token_type> >
{
    lexer() 
      : id("[a-zA-Z0-9]+", 1)
      , st("'[^'\\n]*'", 2)
    {
        self =  id [ 
                    lex::_state = "ST"
                ]
            |   lex::token_def<>(".", 3) [ 
                    lex::_state = "ST" 
                ]
            ;
        
        self("ST") =
                st [ 
                    lex::_state = "INITIAL" 
                ]
            |   lex::token_def<>(".", 4) [(
                    lex::_state = "INITIAL"
                  , lex::_pass = lex::pass_flags::pass_fail 
                )]
            ;
    }
    
    lex::token_def<identifier> id;
    lex::token_def<string_literal> st;
};

typedef lexer::iterator_type token_iterator;

int main()
{
    std::string const s = "foo 'bar'";
    
    content_iterator begin = s.data();
    content_iterator end = s.data() + s.size();
    
    lexer l;
    token_iterator begin2 = l.begin(begin, end, "ST");
    token_iterator end2 = l.end();
    
    std::size_t test_data[] = { 1, 3, 2 };
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
