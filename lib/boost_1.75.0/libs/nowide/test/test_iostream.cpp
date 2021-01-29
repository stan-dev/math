//
//  Copyright (c) 2015 Artyom Beilis (Tonkikh)
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/nowide/iostream.hpp>

#include <boost/nowide/utf/utf.hpp>
#include <limits>
#include <string>

#include "test.hpp"

bool isValidUTF8(const std::string& s)
{
    using namespace boost::nowide::utf;
    for(std::string::const_iterator it = s.begin(); it != s.end();)
    {
        code_point c = utf_traits<char>::decode(it, s.end());
        if(!is_valid_codepoint(c))
            return false;
    }
    return true;
}

void test_main(int argc, char** argv, char**)
{
    const char* example = "Basic letters: \xd7\xa9-\xd0\xbc-\xce\xbd\n"
                          "East Asian Letters: \xe5\x92\x8c\xe5\xb9\xb3\n"
                          "Non-BMP letters: \xf0\x9d\x84\x9e\n"
                          "Invalid UTF-8: `\xFF' `\xd7\xFF' `\xe5\xFF\x8c' `\xf0\x9d\x84\xFF' \n"
                          "\n";

    // If we are using the standard rdbuf we can only put back 1 char
    if(boost::nowide::cin.rdbuf() == std::cin.rdbuf())
    {
        std::cout << "Using std::cin" << std::endl;
        int maxval = 15000;
        for(int i = 0; i < maxval; i++)
        {
            char c = i % 96 + ' ';
            TEST(boost::nowide::cin.putback(c));
            int ci = i % 96 + ' ';
            TEST(boost::nowide::cin.get() == ci);
        }
    } else
    {
        int maxval = 15000;
        for(int i = 0; i < maxval; i++)
        {
            char c = i % 96 + ' ';
            TEST(boost::nowide::cin.putback(c));
        }
        for(int i = maxval - 1; i >= 0; i--)
        {
            int c = i % 96 + ' ';
            TEST(boost::nowide::cin.get() == c);
        }
    }
    boost::nowide::cout << "Normal I/O:" << std::endl;
    boost::nowide::cout << example << std::endl;
    boost::nowide::cerr << example << std::endl;

    boost::nowide::cout << "Flushing each character:" << std::endl;

    for(const char* s = example; *s; s++)
    {
        boost::nowide::cout << *s << std::flush;
        TEST(boost::nowide::cout);
    }

    TEST(boost::nowide::cout);
    TEST(boost::nowide::cerr);
    if(argc == 2 && argv[1] == std::string("-i"))
    {
        boost::nowide::cout << "Input 2 strings" << std::endl;
        std::string v1, v2;
        boost::nowide::cin >> v1 >> v2;
        TEST(boost::nowide::cin);
        TEST(isValidUTF8(v1));
        TEST(isValidUTF8(v2));
        boost::nowide::cout << "First:  " << v1 << std::endl;
        boost::nowide::cout << "Second: " << v2 << std::endl;
        TEST(boost::nowide::cout);

        // Check sync
        boost::nowide::cout << "Input 2 strings\n";
        boost::nowide::cout.flush();
        TEST(boost::nowide::cin >> v1);
        boost::nowide::cin.sync();
        boost::nowide::cout << "First:  " << v1 << std::endl;
        boost::nowide::cout << "2nd string should be ignored. Input 1 more + [ENTER]" << std::endl;
        // And check getline not getting the CR
        TEST(std::getline(boost::nowide::cin, v1));
        TEST(!v1.empty() && v1[v1.size() - 1u] != '\r');
        boost::nowide::cout << "Value:  " << v1 << std::endl;

        boost::nowide::cout << "Press ENTER to exit";
        boost::nowide::cin.clear();
        boost::nowide::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        boost::nowide::cin.get();
    }
}
