# (c) 2019 Qian Wang
# This file is part of nuskybgd released under MIT License.
# See LICENSE file for full details.

"""
Note that when nuskybgd is imported, some checks are performed and information
is stored in nuskybgd.env, including reference to a nuskybgd.caldb.CalDB
object in env._CALDB. If part of this program is run without having first run
this file, the line

env._PASSCHECK = env.check()

Should be run first.
"""

import xspec
# !!Load xspec before astropy!!
# When testing on Mac, using xspec from latest version of HEASOFT (6.26.1) and
# astropy (3.2.1), if astropy is imported first, an error is thrown when
# loading spectrum file, stating that the loaded cfitsio library has SONAME 7.
# HEASOFT 6.26.1 has cfitsio SONAME 8, so probably astropy was compiled with
# cfitsio SONAME 7. At run time the only library in the rpath is HEASOFT's
# latest cfitsio library.

from . import env
import importlib
importlib.import_module('%s.util' % __name__)
importlib.import_module('%s.caldb' % __name__)
importlib.import_module('%s.rmf' % __name__)
importlib.import_module('%s.env' % __name__)

env._PASSCHECK = env.check()
