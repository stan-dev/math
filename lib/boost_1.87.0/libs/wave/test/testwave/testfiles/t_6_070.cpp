/*=============================================================================
    Boost.Wave: A Standard compliant C++ preprocessor library
    http://www.boost.org/

    Copyright (c) 2001-2012 Hartmut Kaiser. Distributed under the Boost
    Software License, Version 1.0. (See accompanying file
    LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

    The tests included in this file were initially taken from the mcpp V2.5
    preprocessor validation suite and were modified to fit into the Boost.Wave
    unit test requirements.
    The original files of the mcpp preprocessor are distributed under the
    license reproduced at the end of this file.
=============================================================================*/

// Tests error reporting: overflow of constant expression in #if directive.

// 14.10:
//E t_6_070.cpp(20): error: integer overflow in preprocessor expression: $E(__TESTWAVE_LONG_MIN__) / - 1
#if __TESTWAVE_LONG_MIN__ / - 1
#endif
