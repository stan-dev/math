from __future__ import absolute_import

import sys
if sys.platform == 'linux2':
    import DLFCN as dl
    flags = sys.getdlopenflags()
    sys.setdlopenflags(dl.RTLD_NOW|dl.RTLD_GLOBAL)
    from . import mpi
    sys.setdlopenflags(flags)
else:
    from . import mpi

