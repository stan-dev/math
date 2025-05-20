//
// Copyright (c) 2022 Evan Lenz (evan@lenzconsulting.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef EXAMPLE_ISSUE_131_HPP
#define EXAMPLE_ISSUE_131_HPP

namespace example {

/** Issue 131

    Properly render multiple variables in a param directive
*/
class issue_131
{
public:
    /** Brief

        @param r,g,b The color components to use
    */
    issue_131(int r, int g, int b);
};

} // example

#endif
