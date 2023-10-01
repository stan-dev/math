//
// Copyright (c) 2022 Alan de Freitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

//[example_route

/*
    This example defines a route for a URL path.
    If the specified route matches and the file
    exists, the example prints its contents to
    standard output.
*/



#include <boost/url/error.hpp>
#include <boost/url/parse.hpp>
#include <boost/url/segments_encoded_ref.hpp>
#include <boost/url/segments_encoded_view.hpp>
#include <boost/url/string_view.hpp>
#include <boost/url/url.hpp>
#include <boost/url/url_view.hpp>
#include <boost/url/static_url.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <iostream>

namespace urls = boost::urls;
namespace fs = boost::filesystem;
using string_view = urls::string_view;

/** Check if a target matches a prefix

    This function checks if the first segments
    of the target match the corresponding prefix
    segments.

    @param target Target segments
    @param prefix Prefix segments
    @return True if target matches prefix
 */
bool match_prefix(
    urls::segments_view target,
    urls::segments_view prefix)
{
    // Trivially reject target that cannot
    // contain the prefix
    if (target.size() < prefix.size())
        return false;

    // Match the prefix segments
    auto it0 = target.begin();
    auto end0 = target.end();
    auto it1 = prefix.begin();
    auto end1 = prefix.end();
    while (
        it0 != end0 &&
        it1 != end1 &&
        *it0 == *it1)
    {
        ++it0;
        ++it1;
    }
    return it1 == end1;
}

/** A static route representing files in a directory

    A route is a URL logical prefix representing
    static files in the specified root directory.

    The `match` function returns the corresponding
    file for a given URL path.
 */
class route
{
public:
    /// Constructor
    route(string_view prefix, fs::path root)
        : prefix_(urls::parse_uri_reference(prefix).value())
        , root_(std::move(root))
    {}

    /// Constructor
    route(urls::url prefix, fs::path root)
        : prefix_(std::move(prefix))
        , root_(std::move(root))
    {}

    /** Match target URL path with a file

        This function attempts to match the target
        URL path with the route prefix.

        If the prefix matches, the target is
        considered to represent a file in the root
        directory. When that happens, the target
        prefix is consumed and other segments are
        appended to the root path.

        The complete file path represented by the
        target is returned as the output parameter
        `result`.

        @param target Target URL path
        @param result An out-parameter holding the resulting path
        @return `true` if target matches the directory
     */
    bool match(
        urls::url_view target,
        fs::path& result)
    {
        if (match_prefix(
                target.segments(),
                static_cast<urls::url_view>(prefix_).segments()))
        {
            result = root_;
            auto segs = target.segments();
            auto it = segs.begin();
            auto end = segs.end();
            std::advance(it, prefix_.segments().size());
            while (it != end)
            {
                auto seg = *it;
                result.append(seg.begin(), seg.end());
                ++it;
            }
            return true;
        }
        return false;
    }

private:
    urls::url prefix_;
    fs::path root_;
};

int
main(int argc, char **argv)
{
    namespace urls = boost::urls;
    namespace fs   = boost::filesystem;

    // Check command line arguments.
    if (argc != 4)
    {
        fs::path exec = argv[0];
        exec = exec.filename();
        std::cerr
            << "Usage: " << exec.c_str()
            << " <target> <prefix> <doc_root>\n"
               "target: path to make a request\n"
               "prefix: url prefix\n"
               "doc_root: dir to look for files\n";
        return EXIT_FAILURE;
    }

    try {
        urls::url target =
            urls::parse_uri_reference(argv[1]).value();
        target.normalize_path();

        std::string prefix = argv[2];
        fs::path root = argv[2];

        if (!fs::is_directory(root))
        {
            std::cerr
                << "Error: " << root
                << " is not a directory\n";
            return EXIT_FAILURE;
        }

        // Create route
        route r(prefix, root);

        // Check if target matches a file
        fs::path result;
        if (r.match(target, result))
        {
            fs::ifstream f(result);
            std::string l;
            while (std::getline(f, l))
                std::cout << l << '\n';
            f.close();
        }
        else
        {
            std::cout
                << "No " << target << " in prefix "
                << prefix << std::endl;
        }
        return EXIT_SUCCESS;
    }
    catch (std::exception &e)
    {
        std::cerr << e.what() << "\n";
        return EXIT_FAILURE;
    }
}

//]
