/*=============================================================================
    Boost.Wave: A Standard compliant C++ preprocessor library
    http://www.boost.org/

    Copyright (c) 2024 Nick Nobles. Distributed under the Boost
    Software License, Version 1.0. (See accompanying file
    LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

// Line directives should be emitted for included files when they begin
// with a single #if/#ifdef/#define directive. Addresses github issue #222.
#include "t_5_040_001.hpp" // #if as first line should emit line directive
#include "t_5_040_002.hpp" // #define as first line should emit line directive
#include "t_5_040_003.hpp" // ensure nested includes emit line directive
#include "t_5_040_005.hpp" // ensure empty initial lines produce line directive

t_5_040_a
#if 1
t_5_040_b
#endif

//R #line 2 "t_5_040_001.hpp"
//R t_5_040_001a
//R #line 6 "t_5_040_001.hpp"
//R t_5_040_001b
//R 
//R t_5_040_001c
//R #line 2 "t_5_040_002.hpp"
//R t_5_040_002
//R #line 2 "t_5_040_004.hpp"
//R t_5_040_004
//R #line 2 "t_5_040_003.hpp"
//R t_5_040_003
//R #line 15 "t_5_040_005.hpp"
//R t_5_040_005
//R #line 17 "t_5_040.cpp"
//R t_5_040_a
//R 
//R t_5_040_b

//H 10: t_5_040.cpp(12): #include "t_5_040_001.hpp"
//H 04: "t_5_040_001.hpp"
//H 05: $S(t_5_040_001.hpp) ($B(t_5_040_001.hpp))
//H 10: t_5_040_001.hpp(1): #if
//H 11: t_5_040_001.hpp(1): #if 1: 1
//H 10: t_5_040_001.hpp(3): #endif
//H 10: t_5_040_001.hpp(7): #if
//H 11: t_5_040_001.hpp(7): #if 1: 1
//H 10: t_5_040_001.hpp(9): #endif
//H 06: 
//H 10: t_5_040.cpp(13): #include "t_5_040_002.hpp"
//H 04: "t_5_040_002.hpp"
//H 05: t_5_040_002.hpp ($B(t_5_040_002.hpp))
//H 10: t_5_040_002.hpp(1): #define
//H 08: t_5_040_002.hpp(1): t_5_040_002_hpp=
//H 06: 
//H 10: t_5_040.cpp(14): #include "t_5_040_003.hpp"
//H 04: "t_5_040_003.hpp"
//H 05: t_5_040_003.hpp ($B(t_5_040_003.hpp))
//H 10: t_5_040_003.hpp(1): #include "t_5_040_004.hpp"
//H 04: "t_5_040_004.hpp"
//H 05: t_5_040_004.hpp ($B(t_5_040_004.hpp))
//H 10: t_5_040_004.hpp(1): #if
//H 11: t_5_040_004.hpp(1): #if 1: 1
//H 10: t_5_040_004.hpp(3): #endif
//H 06: 
//H 06: 
//H 10: t_5_040.cpp(15): #include "t_5_040_005.hpp"
//H 04: "t_5_040_005.hpp"
//H 05: t_5_040_005.hpp ($B(t_5_040_005.hpp))
//H 06: 
//H 10: t_5_040.cpp(18): #if
//H 11: t_5_040.cpp(18): #if 1: 1
//H 10: t_5_040.cpp(20): #endif
