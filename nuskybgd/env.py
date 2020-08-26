# (c) 2019 Qian Wang
# This file is part of nuskybgd released under MIT License.
# See LICENSE file for full details.

import os
from . import env

_PASSCHECK = False
_AUX_ENV = 'NUSKYBGD_AUXIL'
_AUX_DIR = ''
_AUX_FILES = [
    'be.arf',
    'diag.rmf',
    'fcxbA.arf',
    'fcxbB.arf',
    'nomapbgdparams.dat',
    'detA_det1.img.gz',
    'detB_det1.img.gz',
    'ratios.json'
]
_CALDB_ENV = 'CALDB'
_CALDB_CONF_ENV = 'CALDBCONFIG'
_CALDB_PATH = None
_CALDB = None


def check():
    """
    Check for configuration problems.

    - ``NUSKYBGD_AUXIL`` shell environment variable is set. The name of this
      variable is defined by ``nuskybgd.env._AUX_ENV``.
    - Auxiliary files listed in ``nuskybgd.env._AUX_FILES`` are present.
    - ``CALDB`` shell environment variable is set. The name of this variable
      is defined by ``nuskybgd.env._CALDB_ENV``.
    - ``CALDBCONFIG`` shell environment variable is set. The name of this
      variable is defined by ``nuskybgd.env._CALDB_CONF_ENV``.
    - An instance of ``nuskybgd.caldb.CalDB()`` is created successfully using
      the path in ``CALDB`` environment variable.

    If these conditions are met, ``nuskybgd.env._PASSCHECK`` is set to
    ``True`` and subsequent calls to ``nuskybgd.env.block()`` return ``False``.
    """
    from . import caldb

    fail = False

    # Auxil directory
    if env._AUX_ENV not in os.environ:
        print('Please set the %s environment variable first.' % env._AUX_ENV)
        fail = True
    else:
        env._AUX_DIR = os.environ[env._AUX_ENV]

        # Check auxil dir is OK...
        for _ in env._AUX_FILES:
            if not os.path.exists('%s/%s' % (env._AUX_DIR, _)):
                print('Error: %s not in auxil directory!' % _)
                fail = True

    # CALDB check
    if env._CALDB_ENV not in os.environ:
        print('Please set the %s environment variable first.' % env._CALDB_ENV)
        fail = True
    elif env._CALDB_CONF_ENV not in os.environ:
        print('Please set the %s environment variable first.' % env._CALDB_CONF_ENV)
        fail = True
    else:
        env._CALDB_PATH = os.environ[env._CALDB_ENV]
        cal = caldb.CalDB(env._CALDB_PATH)
        if cal._index is None:
            fail = True

    if fail is True:
        print('Initial check did not pass!')
    else:
        env._CALDB = cal

    return not fail


def block():
    """
    Return ``True`` if ``env.check()`` did not pass, and prints an info message.
    """
    if env._PASSCHECK is False:
        print('Errors encountered when checking setup for nuskybgd.')
        print('Please fix the problems before continuing.')
        return True
    else:
        return False
