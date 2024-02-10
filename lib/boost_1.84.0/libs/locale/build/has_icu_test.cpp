// Copyright (c) 2010 John Maddock
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <unicode/coll.h>
#include <unicode/locid.h>
#include <unicode/uchar.h>
#include <unicode/utypes.h>
#include <unicode/uversion.h>

int main()
{
    icu::Locale loc;
    UErrorCode err = U_ZERO_ERROR;
    UChar32 c = ::u_charFromName(U_UNICODE_CHAR_NAME, "GREEK SMALL LETTER ALPHA", &err);
    return err;
}
