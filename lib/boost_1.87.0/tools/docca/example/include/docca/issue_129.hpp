//
// Copyright (c) 2022 Evan Lenz (evan@lenzconsulting.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef EXAMPLE_ISSUE_129_HPP
#define EXAMPLE_ISSUE_129_HPP

namespace example {

/** Issue 129

    Render replacement text in code in comments
*/
class issue_129
{
public:
    /** Brief

        @par Value Type
        @code
        using value_type = __see_below__;
        using value_type = __implementation_defined__;
        @endcode
    */
    issue_129(int foo);
};

} // example

#endif
