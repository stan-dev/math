#!/usr/bin/env python
from __future__ import print_function
import sys
import re

"""
Removes macros that might otherwise cause confusion for CMock.
"""

if len(sys.argv) != 3:
    print("Usage: strip_defines.py <input> <output>")
    sys.exit(1)
with open(sys.argv[1], 'r') as inf:
    text = inf.read()
text = re.sub(r'\bCL_(?:API_ENTRY|API_SUFFIX|EXT|CALLBACK)[A-Z0-9_]*\b', '', text)
text2 = re.sub(r'\*\[\]', '\*\*', text)
with open(sys.argv[2], 'w') as outf:
    outf.write(text2)
