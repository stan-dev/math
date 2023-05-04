//
// Copyright (c) 2022 alandefreitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//

//[example_suffix_list

/*
    The Public Suffix List is intended to
    enumerate all domain suffixes controlled
    by registrars.
    https://publicsuffix.org/list/
*/

#include <boost/url/url.hpp>
#include <boost/url/parse.hpp>
#include <fstream>
#include <iostream>

namespace urls = boost::urls;

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cout << argv[0] << "\n";
        std::cout << "Usage: suffix_list <url> <suffix_list>\n"
                     "options:\n"
                     "    <url>:              A valid url (required)\n"
                     "    <suffix_list>:      File with the public suffix list (default: public_suffix_list.dat)\n"
                     "examples:\n"
                     "suffix_list \"www.example.com\" \"public_suffix_list.dat\"\n";
        return EXIT_FAILURE;
    }

    urls::url u = urls::parse_uri(argv[1]).value();

    std::string filename =
        argc < 3 ?
            "public_suffix_list.dat" :
            argv[2];
    std::ifstream fin(filename);
    if (fin.bad())
    {
        std::cerr << "Cannot open " << filename << "\n";
        return EXIT_FAILURE;
    }

    std::string comment;
    std::string suffix;
    std::string line;
    while (std::getline(fin, line))
    {
        urls::string_view sv(line);
        if (sv.starts_with("//") && suffix.empty())
        {
            comment.append(sv.substr(3));
            comment.append("\n");
        }
        else if (line.empty())
        {
            if (suffix.empty())
                comment.clear();
            else
                break;
        }
        else if (u.encoded_host().ends_with(line))
        {
            suffix = line;
        }
    }

    std::cout <<
        "url:    \n" << u                << "\n\n"
        "host:   \n" << u.encoded_host() << "\n\n"
        "suffix: \n" << suffix           << "\n\n"
        "detail: \n" << comment          << "\n\n";

    return EXIT_SUCCESS;
}

//]
