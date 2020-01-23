#ifndef B2_CONFIG_H
#define B2_CONFIG_H

/*
Copyright 2002-2018 Rene Rivera.
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or copy at
http://www.boost.org/LICENSE_1_0.txt)
*/

#define OPT_HEADER_CACHE_EXT 1
#define OPT_GRAPH_DEBUG_EXT 1
#define OPT_SEMAPHORE 1
#define OPT_AT_FILES 1
#define OPT_DEBUG_PROFILE 1
#define JAM_DEBUGGER 1
#define OPT_FIX_TARGET_VARIABLES_EXT 1
#define OPT_IMPROVED_PATIENCE_EXT 1

// Autodetect various operating systems..

#if defined(_WIN32) || defined(_WIN64) || \
    defined(__WIN32__) || defined(__TOS_WIN__) || \
    defined(__WINDOWS__)
    #define NT 1
#endif

#if defined(__VMS) || defined(__VMS_VER)
    #if !defined(VMS)
        #define VMS 1
    #endif
#endif

#endif
