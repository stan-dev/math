/*
Copyright 2020 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#ifndef B2_STARTUP_H
#define B2_STARTUP_H

#include "config.h"
#include "frames.h"

namespace b2
{
    namespace startup
    {
        void load_builtins();
        LIST *builtin_boost_build(FRAME *frame, int flags);
        bool bootstrap(FRAME *frame);
    }
} // namespace b2

#endif
