//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/CPPAllinace/url
//

Correct hyphenation:

Nouns: never hyphenated
1. The percent encoding
2. Each percent escape

Adjectives: hyphenated when the noun comes after
3. Returns a percent-escaped string
4. Removes percent-encoding from the string

Adjectives: not hyphenated when the noun comes before
5. The string `s` must be percent encoded.
6. The string may not be percent escaped.



// https://eel.is/c++draft/structure.specifications#3

    /** {brief}

        {description}

        Any percent-escapes in the string are
        decoded first.

        The returned string may contain
        percent escapes.

        Reserved characters in the string are
        percent-escaped in the result.

        Escapes in the string are preserved,
        and reserved characters in the string
        are percent-escaped in the result.

        The comparison is performed as if all
        escaped characters were decoded first.

        @note
        The interpretation of userinfo as
        individual user and password components
        is scheme-dependent. Transmitting
        passwords in URLs is deprecated.

        @par Example

        @par Constraints    (SFINAE)
        @par Mandates       (static_assert)
        @par Preconditions  (assert)
        @par Effects
        @par Synchronization
        @par Postconditions
        @par Complexity

        @par Exception Safety
        Strong guarantee.
        Basic guarantee.
        Calls to allocate may throw.
        Exceptions thrown on invalid input.

        @throw system_error
        `s` contains an invalid percent-encoding.

        @tparam
        @return
        @param

        @par BNF
        @par Specification
        @li

        @see
            @ref
    */

//------------------------------------------------

/** {brief}

    {description}
    <!-- Assume the user knows C++ properly -->
    <!-- Assume the user has read the whole exposition -->
    <!-- Avoid too much symbolism -->
    <!-- Avoid repeating words -->
    <!-- Use the present tense -->
    <!-- Lists should use noun-phrases -->

    @par Example
    @code
    {example}
    // make it short
    @endcode

    @par Mandates

    @par Preconditions

    @par Effects

    @note {text}

    @par Exception Safety
    <!-- exception guarantee + explanation -->

    <!-- Exception guarantee -->
    Throws nothing. <!-- the function never throws exceptions -->
    Strong guarantee. <!-- If the function throws an exception, the state of the program is rolled back -->
    Basic guarantee. <!-- If the function throws an exception, the program is in a valid state -->
    No guarantee. <!-- If the function throws an exception, the program may not be in a valid state -->

    <!-- Use explanation stock phrases if applicable: -->
    <!-- Calls to allocate may throw. -->
    <!-- Exceptions thrown on invalid input. -->
    <!-- Exceptions thrown on excessive input length. -->

    @tparam {name} {description}.

    @return {description}.

    @param {name} {description}.

    <!-- Use stock phrases when applicable -->
    <!-- @param begin An iterator to the beginning of the range -->
    <!-- @param end An iterator to the end of the range -->
    <!-- @param ec Set to the error, if any occurred -->

    @throws {exception} {condition}.

    <!-- Use stock phrases when applicable -->
    <!-- @throws boost::system::system_error Thrown on failure. -->

    @par Thread Safety
    {description}

    @par Specification
    @li <a href="">text (rfc#)</a>

    @par References
    @li <a href="">text</a>

    @see
        @ref {refid},
        @ref {refid},
        @ref {refid}.

*/