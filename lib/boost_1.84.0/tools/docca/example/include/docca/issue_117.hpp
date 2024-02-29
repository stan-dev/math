//
// Copyright (c) 2022 Evan Lenz (evan@lenzconsulting.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef EXAMPLE_ISSUE_117_HPP
#define EXAMPLE_ISSUE_117_HPP

namespace example {

/** Issue 117

    Properly detect (and omit) empty description sections, even when
    Doxygen puts the param list inside the detaileddescription.
*/
class issue_117
{
public:
    /** Construct from an array of bytes.

        @param bytes The value to construct from.
    */
    issue_117(
        bytes_type const& bytes) noexcept;
};

} // example

#endif
