//
// Copyright (c) 2009-2011 Artyom Beilis (Tonkikh)
// Copyright (c) 2023 Alexander Grund
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

//
// ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
//
// BIG FAT WARNING FOR Microsoft Visual Studio Users
//
// YOU NEED TO CONVERT THIS SOURCE FILE ENCODING TO UTF-8 WITH BOM ENCODING.
//
// Unfortunately MSVC understands that the source code is encoded as
// UTF-8 only if you add useless BOM in the beginning.
//
// So, before you compile "wide" examples with MSVC, please convert them to text
// files with BOM. There are two very simple ways to do it:
//
// 1. Open file with Notepad and save it from there. It would convert
//    it to file with BOM.
// 2. In Visual Studio go File->Advances Save Options... and select
//    Unicode (UTF-8  with signature) Codepage 65001
//
// Note: once converted to UTF-8 with BOM, this source code would not
// compile with other compilers, because no-one uses BOM with UTF-8 today
// because it is absolutely meaningless in context of UTF-8.
//
// ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
//
#include <boost/locale.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <iostream>

int main()
{
    if(boost::locale::localization_backend_manager::global().get_all_backends().at(0) != "icu")
        std::wcout << L"Need ICU support for this example!\nConversion below will likely be wrong!\n";

    // Create system default locale
    boost::locale::generator gen;
    std::locale loc = gen("");
    std::locale::global(loc);
    std::wcout.imbue(loc);

    // This is needed to prevent the C stdio library from
    // converting strings to narrow on some platforms
    std::ios_base::sync_with_stdio(false);

    std::wcout << L"Correct case conversion can't be done by simple, character by character conversion\n";
    std::wcout << L"because case conversion is context sensitive and not a 1-to-1 conversion.\n";
    std::wcout << L"For example:\n";
    const std::wstring gruessen(L"grüßen");
    std::wcout << L"   German " << gruessen << " would be incorrectly converted to " << boost::to_upper_copy(gruessen);
    std::wcout << ", while Boost.Locale converts it to " << boost::locale::to_upper(gruessen) << std::endl
               << L"     where ß is replaced with SS.\n";
    const std::wstring greek(L"ὈΔΥΣΣΕΎΣ");
    std::wcout << L"   Greek " << greek << " would be incorrectly converted to " << boost::to_lower_copy(greek);
    std::wcout << ", while Boost.Locale correctly converts it to " << boost::locale::to_lower(greek) << std::endl
               << L"     where Σ is converted to σ or to ς, according to position in the word.\n";
    std::wcout
      << L"Such type of conversion just can't be done using std::toupper/boost::to_upper* that work on character "
         L"by character base.\n"
         L"Also std::toupper is not fully applicable when working with variable character length like UTF-8 or UTF-16\n"
         L"limiting the correct behavior to BMP or ASCII only\n";
}

// boostinspect:noascii
