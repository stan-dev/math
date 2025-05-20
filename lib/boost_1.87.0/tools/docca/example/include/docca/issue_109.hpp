//
// Copyright (c) 2022 Evan Lenz (evan@lenzconsulting.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef EXAMPLE_ISSUE_109_HPP
#define EXAMPLE_ISSUE_109_HPP

namespace example {

/** Issue 109

    Properly detect (and omit) empty description sections when
    Doxygen puts simplesect elements inside the detaileddescription.
*/
class issue_109
{
public:
    /** Constructor

        @par Test section

            Test paragraph
    */
    issue_109();
};

} // example

#endif
