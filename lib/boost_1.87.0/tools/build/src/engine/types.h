/*  Copyright 2019-2023 René Ferdinand Rivera Morell
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE.txt or
 * https://www.bfgroup.xyz/b2/LICENSE.txt)
 */

#ifndef B2_TYPES_H
#define B2_TYPES_H

#include "config.h"

#include <string>

#if B2_USE_STD_THREADS
#include <mutex>
#endif

namespace b2 {
using string_t = std::string;
using int_t = int;
using uint_t = unsigned int;

#if B2_USE_STD_THREADS

using mutex_t = std::mutex;
using scope_lock_t = std::unique_lock<std::mutex>;

#else

struct mutex_t
{};
struct scope_lock_t
{
	inline scope_lock_t(mutex_t &) {}
};

#endif

} // namespace b2

#endif
