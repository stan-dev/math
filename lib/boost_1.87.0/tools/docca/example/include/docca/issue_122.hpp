//
// Copyright (c) 2022 Evan Lenz (evan@lenzconsulting.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef EXAMPLE_ISSUE_122_HPP
#define EXAMPLE_ISSUE_122_HPP

namespace example {

/** Issue 122

    Properly detect and render deleted members.
*/
class issue_122
{
public:
    /** Construct something.

        @param foo The value to construct from.
    */
    issue_122(int foo);
    issue_122(float foo) = delete;
};

} // example

#endif
