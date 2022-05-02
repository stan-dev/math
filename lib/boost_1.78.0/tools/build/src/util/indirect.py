# Status: minimally ported. This module is not supposed to be used much
# with Boost.Build/Python.
#
# Copyright 2003 Dave Abrahams
# Copyright 2003 Vladimir Prus
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)

from b2.util import call_jam_function, bjam_signature

def call(*args):
    a1 = args[0]
    name = a1[0]
    a1tail = a1[1:]
    call_jam_function(name, *((a1tail,) + args[1:]))
