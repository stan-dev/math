#ifndef B2_CONFIG_H
#define B2_CONFIG_H

/*
Copyright 2002-2023 Rene Rivera.
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or copy at
https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

// Sigh, Widows.
#define NOMINMAX

#define OPT_HEADER_CACHE_EXT 1
#define OPT_GRAPH_DEBUG_EXT 1
#define OPT_SEMAPHORE 1
#define OPT_AT_FILES 1
#define OPT_DEBUG_PROFILE 1
#define JAM_DEBUGGER 1
#define OPT_FIX_TARGET_VARIABLES_EXT 1
#define OPT_IMPROVED_PATIENCE_EXT 1
// #define BJAM_NO_MEM_CACHE 1

// Autodetect various operating systems..

#if defined(_WIN32) || defined(_WIN64) || defined(__WIN32__) \
	|| defined(__TOS_WIN__) || defined(__WINDOWS__)
#define NT 1
#endif

#if defined(__VMS) || defined(__VMS_VER)
#if !defined(VMS)
#define VMS 1
#endif
#endif

// To work around QEMU failures on mixed mode situations (32 vs 64) we need to
// enable partial LFS support in system headers. And we need to do this before
// any system headers are included.

#if !defined(NT) && !defined(VMS)
#define _FILE_OFFSET_BITS 64
#endif

// Account for incomplete C++ standard implementations.

#ifndef B2_NOEXCEPT
#define B2_NOEXCEPT noexcept
#endif

// Indicate if we can use std::thread and friends.
#ifndef B2_USE_STD_THREADS
#define B2_USE_STD_THREADS 1
#endif

#include <cassert>

namespace b2 {

inline bool ensure(bool v)
{
#if defined(B2_DEBUG) && B2_DEBUG
	assert(v);
#endif
	return v;
}

template <typename T>
T * ensure_valid(T * v)
{
#if defined(B2_DEBUG) && B2_DEBUG
	assert(v != nullptr);
#endif
	return v;
}

} // namespace b2

#endif
