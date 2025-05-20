//
// Copyright (c) 2022 Evan Lenz (evan@lenzconsulting.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef EXAMPLE_ISSUE_133_HPP
#define EXAMPLE_ISSUE_133_HPP

namespace example {

/** Issue 133

    Render tables supported by Doxygen
*/
class issue_133
{
public:
    /** Table example

    First Header  | Second Header
    ------------- | -------------
    Content Cell  | Content Cell
    Content Cell  | Content Cell
     
    | Right | Center | Left  |
    | ----: | :----: | :---- |
    | 10    | 10     | 10    |
    | 1000  | 1000   | 1000  |

    | Right | Center | Left  |
    | ----: | :----: | :---- |
    | 10    | 10     | 10    |
    | ^     | 1000   | 1000  |

    | Right | Center | Left  |
    | ----: | :----: | :---- |
    | 10    | 10     | 10    |
    | 1000  |||
    
    */
    issue_133();
};

} // example

#endif
