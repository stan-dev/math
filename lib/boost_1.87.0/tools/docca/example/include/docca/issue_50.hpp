//
// Copyright (c) 2016-2019 Vinnie Falco (vinnie dot falco at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/beast
//

#ifndef EXAMPLE_ISSUE_50
#define EXAMPLE_ISSUE_50

namespace example {

/** A basic issue.

    Lorum ipsum...
*/
template<class Allocator>
class basic_issue_50
{
};

/** A typical issue 50.
*/
using issue_50 = basic_issue_50<std::allocator<char>>;

} // example

#endif
