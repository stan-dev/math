//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
// Copyright (c) 2023 Alexander Grund
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/locale.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <iostream>

int main()
{
    if(boost::locale::localization_backend_manager::global().get_all_backends().at(0) != "icu")
        std::cout << "Need ICU support for this example!\nConversion below will likely be wrong!\n";

    // Create system default locale
    boost::locale::generator gen;
    std::locale loc = gen("");
    std::locale::global(loc);
    std::cout.imbue(loc);

    // This is needed to prevent the C stdio library from
    // converting strings to narrow on some platforms
    std::ios_base::sync_with_stdio(false);

    std::cout << "Correct case conversion can't be done by simple, character by character conversion\n";
    std::cout << "because case conversion is context sensitive and not a 1-to-1 conversion.\n";
    std::cout << "For example:\n";
    const std::string gruessen("grüßen");
    std::cout << "   German " << gruessen << " would be incorrectly converted to " << boost::to_upper_copy(gruessen);
    std::cout << ", while Boost.Locale converts it to " << boost::locale::to_upper(gruessen) << std::endl
              << "     where ß is replaced with SS.\n";
    const std::string greek("ὈΔΥΣΣΕΎΣ");
    std::cout << "   Greek " << greek << " would be incorrectly converted to " << boost::to_lower_copy(greek);
    std::cout << ", while Boost.Locale correctly converts it to " << boost::locale::to_lower(greek) << std::endl
              << "     where Σ is converted to σ or to ς, according to position in the word.\n";
    std::cout
      << "Such type of conversion just can't be done using std::toupper/boost::to_upper* that work on character "
         "by character base.\n"
         "Also std::toupper is not fully applicable when working with variable character length like UTF-8 or UTF-16\n"
         "limiting the correct behavior to BMP or ASCII only\n";
}

// boostinspect:noascii
