//
//  Copyright (c) 2015 Artyom Beilis (Tonkikh)
//  Copyright (c) 2019 Alexander Grund
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/nowide/fstream.hpp>

#include <boost/nowide/convert.hpp>
#include <boost/nowide/cstdio.hpp>
#include <fstream>
#include <iostream>
#include <string>

#include "test.hpp"

namespace nw = boost::nowide;

void make_empty_file(const char* filepath)
{
    nw::ofstream f(filepath, std::ios_base::out | std::ios::trunc);
    TEST(f);
}

bool file_exists(const char* filepath)
{
    FILE* f = nw::fopen(filepath, "r");
    if(f)
    {
        std::fclose(f);
        return true;
    } else
        return false;
}

std::string read_file(const char* filepath, bool binary_mode = false)
{
    FILE* f = nw::fopen(filepath, binary_mode ? "rb" : "r");
    TEST(f);
    std::string content;
    int c;
    while((c = std::fgetc(f)) != EOF)
        content.push_back(static_cast<char>(c));
    std::fclose(f);
    return content;
}

void test_with_different_buffer_sizes(const char* filepath)
{
    /* Important part of the standard for mixing input with output:
       However, output shall not be directly followed by input without an intervening call to the fflush function
       or to a file positioning function (fseek, fsetpos, or rewind),
       and input shall not be directly followed by output without an intervening call to a file positioning function,
       unless the input operation encounters end-of-file.
    */
    for(int i = -1; i < 16; i++)
    {
        std::cout << "Buffer size = " << i << std::endl;
        char buf[16];
        nw::fstream f;
        // Different conditions when setbuf might be called: Usually before opening a file is OK
        if(i >= 0)
            f.rdbuf()->pubsetbuf((i == 0) ? NULL : buf, i);
        f.open(filepath, std::ios::in | std::ios::out | std::ios::trunc | std::ios::binary);
        TEST(f);
        // Add 'abcdefg'
        TEST(f.put('a'));
        TEST(f.put('b'));
        TEST(f.put('c'));
        TEST(f.put('d'));
        TEST(f.put('e'));
        TEST(f.put('f'));
        TEST(f.put('g'));
        // Read first char
        TEST(f.seekg(0));
        TEST(f.get() == 'a');
        TEST(f.gcount() == 1u);
        // Skip next char
        TEST(f.seekg(1, std::ios::cur));
        TEST(f.get() == 'c');
        TEST(f.gcount() == 1u);
        // Go back 1 char
        TEST(f.seekg(-1, std::ios::cur));
        TEST(f.get() == 'c');
        TEST(f.gcount() == 1u);

        // Test switching between read->write->read
        // case 1) overwrite, flush, read
        TEST(f.seekg(1));
        TEST(f.put('B'));
        TEST(f.flush()); // Flush when changing out->in
        TEST(f.get() == 'c');
        TEST(f.gcount() == 1u);
        TEST(f.seekg(1));
        TEST(f.get() == 'B');
        TEST(f.gcount() == 1u);
        // case 2) overwrite, seek, read
        TEST(f.seekg(2));
        TEST(f.put('C'));
        TEST(f.seekg(3)); // Seek when changing out->in
        TEST(f.get() == 'd');
        TEST(f.gcount() == 1u);

        // Check that sequence from start equals expected
        TEST(f.seekg(0));
        TEST(f.get() == 'a');
        TEST(f.get() == 'B');
        TEST(f.get() == 'C');
        TEST(f.get() == 'd');
        TEST(f.get() == 'e');

        // Putback after flush is implementation defined
        // Boost.Nowide: Works
#if BOOST_NOWIDE_USE_FILEBUF_REPLACEMENT
        TEST(f << std::flush);
        TEST(f.putback('e'));
        TEST(f.putback('d'));
        TEST(f.get() == 'd');
        TEST(f.get() == 'e');
#endif
        // Rest of sequence
        TEST(f.get() == 'f');
        TEST(f.get() == 'g');
        TEST(f.get() == EOF);

        // Put back until front of file is reached
        f.clear();
        TEST(f.seekg(1));
        TEST(f.get() == 'B');
        TEST(f.putback('B'));
        // Putting back multiple chars is not possible on all implementations after a seek/flush
#if BOOST_NOWIDE_USE_FILEBUF_REPLACEMENT
        TEST(f.putback('a'));
        TEST(!f.putback('x')); // At beginning of file -> No putback possible
        // Get characters that were putback to avoid MSVC bug https://github.com/microsoft/STL/issues/342
        f.clear();
        TEST(f.get() == 'a');
#endif
        TEST(f.get() == 'B');
        f.close();
        TEST(nw::remove(filepath) == 0);
    }
}

void test_close(const char* filepath)
{
    const std::string filepath2 = std::string(filepath) + ".2";
    // Make sure file does not exist yet
    TEST(!file_exists(filepath2.c_str()) || nw::remove(filepath2.c_str()) == 0);
    TEST(!file_exists(filepath2.c_str()));
    nw::filebuf buf;
    TEST(buf.open(filepath, std::ios_base::out) == &buf);
    TEST(buf.is_open());
    // Opening when already open fails
    TEST(buf.open(filepath2.c_str(), std::ios_base::out) == NULL);
    // Still open
    TEST(buf.is_open());
    TEST(buf.close() == &buf);
    // Failed opening did not create file
    TEST(!file_exists(filepath2.c_str()));
    // But it should work now:
    TEST(buf.open(filepath2.c_str(), std::ios_base::out) == &buf);
    TEST(buf.close() == &buf);
    TEST(file_exists(filepath2.c_str()));
    TEST(nw::remove(filepath) == 0);
    TEST(nw::remove(filepath2.c_str()) == 0);
}

template<typename IFStream, typename OFStream>
void test_flush(const char* filepath)
{
    OFStream fo(filepath, std::ios_base::out | std::ios::trunc);
    TEST(fo);
    std::string curValue;
    for(int repeat = 0; repeat < 2; repeat++)
    {
        for(size_t len = 1; len <= 1024; len *= 2)
        {
            char c = static_cast<char>(len % 13 + repeat + 'a'); // semi-random char
            std::string input(len, c);
            fo << input;
            curValue += input;
            TEST(fo.flush());
            std::string s;
            // Note: Flush on read area is implementation defined, so check whole file instead
            IFStream fi(filepath);
            TEST(fi >> s);
            // coverity[tainted_data]
            TEST(s == curValue);
        }
    }
}

void test_ofstream_creates_file(const char* filename)
{
    TEST(!file_exists(filename) || nw::remove(filename) == 0);
    TEST(!file_exists(filename));
    // Ctor
    {
        nw::ofstream fo(filename);
        TEST(fo);
    }
    TEST(file_exists(filename));
    TEST(read_file(filename).empty());
    TEST(nw::remove(filename) == 0);
    // Open
    {
        nw::ofstream fo;
        fo.open(filename);
        TEST(fo);
    }
    TEST(file_exists(filename));
    TEST(read_file(filename).empty());
    TEST(nw::remove(filename) == 0);
}

// Create filename file with content "test\n"
void test_ofstream_write(const char* filename)
{
    // char* ctor
    {
        nw::ofstream fo(filename);
        TEST(fo << "test" << 2 << std::endl);
    }
    // char* open
    TEST(read_file(filename) == "test2\n");
    TEST(nw::remove(filename) == 0);
    {
        nw::ofstream fo;
        fo.open(filename);
        TEST(fo << "test" << 2 << std::endl);
    }
    TEST(read_file(filename) == "test2\n");
    TEST(nw::remove(filename) == 0);
    // string ctor
    {
        std::string name = filename;
        nw::ofstream fo(name);
        TEST(fo << "test" << 2 << std::endl);
    }
    TEST(read_file(filename) == "test2\n");
    TEST(nw::remove(filename) == 0);
    // string open
    {
        nw::ofstream fo;
        fo.open(std::string(filename));
        TEST(fo << "test" << 2 << std::endl);
    }
    TEST(read_file(filename) == "test2\n");
    TEST(nw::remove(filename) == 0);
    // Binary mode
    {
        nw::ofstream fo(filename, std::ios::binary);
        TEST(fo << "test" << 2 << std::endl);
    }
    TEST(read_file(filename, true) == "test2\n");
    TEST(nw::remove(filename) == 0);
    // At end
    {
        {
            nw::ofstream fo(filename);
            TEST(fo << "test" << 2 << std::endl);
        }
        nw::ofstream fo(filename, std::ios::ate | std::ios::in);
        fo << "second" << 2 << std::endl;
    }
    TEST(read_file(filename) == "test2\nsecond2\n");
    TEST(nw::remove(filename) == 0);
}

void test_ifstream_open_read(const char* filename)
{
    // Create test file
    {
        nw::ofstream fo(filename);
        TEST(fo << "test" << std::endl);
    }

    // char* Ctor
    {
        nw::ifstream fi(filename);
        TEST(fi);
        std::string tmp;
        TEST(fi >> tmp);
        TEST(tmp == "test");
    }
    // char* open
    {
        nw::ifstream fi;
        fi.open(filename);
        TEST(fi);
        std::string tmp;
        TEST(fi >> tmp);
        TEST(tmp == "test");
    }
    // string ctor
    {
        std::string name = filename;
        nw::ifstream fi(name);
        TEST(fi);
        std::string tmp;
        TEST(fi >> tmp);
        TEST(tmp == "test");
    }
    // string open
    {
        nw::ifstream fi;
        fi.open(std::string(filename));
        TEST(fi);
        std::string tmp;
        TEST(fi >> tmp);
        TEST(tmp == "test");
    }
    // Binary mode
    {
        nw::ifstream fi(filename, std::ios::binary);
        TEST(fi);
        std::string tmp;
        TEST(fi >> tmp);
        TEST(tmp == "test");
    }
    // At end
    {
        // Need binary file or position check might be throw off by newline conversion
        {
            nw::ofstream fo(filename, nw::fstream::binary);
            TEST(fo << "test");
        }
        nw::ifstream fi(filename, nw::fstream::ate | nw::fstream::binary);
        TEST(fi);
        TEST(fi.tellg() == std::streampos(4));
        fi.seekg(-2, std::ios_base::cur);
        std::string tmp;
        TEST(fi >> tmp);
        TEST(tmp == "st");
    }
    // Fail on non-existing file
    TEST(nw::remove(filename) == 0);
    {
        nw::ifstream fi(filename);
        TEST(!fi);
    }
}

void test_fstream(const char* filename)
{
    const std::string sFilename = filename;
    TEST(!file_exists(filename) || nw::remove(filename) == 0);
    TEST(!file_exists(filename));
    // Fail on non-existing file
    {
        nw::fstream f(filename);
        TEST(!f);
        nw::fstream f2(sFilename);
        TEST(!f2);
    }
    {
        nw::fstream f;
        f.open(filename);
        TEST(!f);
        f.open(sFilename);
        TEST(!f);
    }
    TEST(!file_exists(filename));
    // Create empty file (Ctor)
    {
        nw::fstream f(filename, std::ios::out);
        TEST(f);
    }
    TEST(read_file(filename).empty());
    // Char* ctor
    {
        nw::fstream f(filename);
        TEST(f);
        TEST(f << "test");
        std::string tmp;
        TEST(f.seekg(0));
        TEST(f >> tmp);
        TEST(tmp == "test");
    }
    TEST(read_file(filename) == "test");
    // String ctor
    {
        nw::fstream f(sFilename);
        TEST(f);
        TEST(f << "string_ctor");
        std::string tmp;
        TEST(f.seekg(0));
        TEST(f >> tmp);
        TEST(tmp == "string_ctor");
    }
    TEST(read_file(filename) == "string_ctor");
    TEST(nw::remove(filename) == 0);
    // Create empty file (open)
    {
        nw::fstream f;
        f.open(filename, std::ios::out);
        TEST(f);
    }
    TEST(read_file(filename).empty());
    // Open
    {
        nw::fstream f;
        f.open(filename);
        TEST(f);
        TEST(f << "test");
        std::string tmp;
        TEST(f.seekg(0));
        TEST(f >> tmp);
        TEST(tmp == "test");
    }
    TEST(read_file(filename) == "test");
    // Ctor existing file
    {
        nw::fstream f(filename);
        TEST(f);
        std::string tmp;
        TEST(f >> tmp);
        TEST(tmp == "test");
        TEST(f.eof());
        f.clear();
        TEST(f << "second");
    }
    TEST(read_file(filename) == "testsecond");
    // Trunc & binary
    {
        nw::fstream f(filename, std::ios::in | std::ios::out | std::ios::trunc | std::ios::binary);
        TEST(f);
        TEST(f << "test2");
        std::string tmp;
        TEST(f.seekg(0));
        TEST(f >> tmp);
        TEST(tmp == "test2");
    }
    TEST(read_file(filename) == "test2");
    // Reading in write mode fails (existing file!)
    {
        nw::fstream f(filename, std::ios::out);
        std::string tmp;
        TEST(!(f >> tmp));
        f.clear();
        TEST(f << "foo");
        TEST(f.seekg(0));
        TEST(!(f >> tmp));
    }
    TEST(read_file(filename) == "foo");
    // Writing in read mode fails (existing file!)
    {
        nw::fstream f(filename, std::ios::in);
        TEST(!(f << "bar"));
        f.clear();
        std::string tmp;
        TEST(f >> tmp);
        TEST(tmp == "foo");
    }
    TEST(read_file(filename) == "foo");
    TEST(nw::remove(filename) == 0);
}

template<typename T>
bool is_open(T& stream)
{
    // There are const and non const versions of is_open, so test both
    TEST(stream.is_open() == const_cast<const T&>(stream).is_open());
    return stream.is_open();
}

template<typename T>
void do_test_is_open(const char* filename)
{
    T f;
    TEST(!is_open(f));
    f.open(filename);
    TEST(f);
    TEST(is_open(f));
    f.close();
    TEST(f);
    TEST(!is_open(f));
}

void test_is_open(const char* filename)
{
    // Note the order: Output before input so file exists
    do_test_is_open<nw::ofstream>(filename);
    do_test_is_open<nw::ifstream>(filename);
    do_test_is_open<nw::fstream>(filename);
    TEST(nw::remove(filename) == 0);
}

void test_main(int, char** argv, char**)
{
    const std::string exampleFilename = std::string(argv[0]) + "-\xd7\xa9-\xd0\xbc-\xce\xbd.txt";
    std::cout << "Testing fstream" << std::endl;
    test_ofstream_creates_file(exampleFilename.c_str());
    test_ofstream_write(exampleFilename.c_str());
    test_ifstream_open_read(exampleFilename.c_str());
    test_fstream(exampleFilename.c_str());
    test_is_open(exampleFilename.c_str());

    std::cout << "Complex IO" << std::endl;
    test_with_different_buffer_sizes(exampleFilename.c_str());

    std::cout << "filebuf::close" << std::endl;
    test_close(exampleFilename.c_str());

    std::cout << "Flush - Sanity Check" << std::endl;
    test_flush<std::ifstream, std::ofstream>(exampleFilename.c_str());
    std::cout << "Flush - Test" << std::endl;
    test_flush<nw::ifstream, nw::ofstream>(exampleFilename.c_str());
}
