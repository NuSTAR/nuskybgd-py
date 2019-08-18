import os

_PASSCHECK = False
_AUX_ENV = 'NUSKYBGD_AUXIL'
_AUX_DIR = ''
_AUX_FILES = [
    'be.arf',
    'diag.rmf',
    'fcxbA.arf',
    'fcxbB.arf',
    'nomapbgdparams.dat',
    'ratios.json'
]
_CALDB_ENV = 'CALDB'
_CALDB_CONF_ENV = 'CALDBCONFIG'


def check():
    """
    Check for configuration problems.
    """
    from . import caldb

    fail = False

    # Auxil directory
    if _AUX_ENV not in os.environ:
        print('Please set the %s environment variable first.' % _AUX_ENV)
        fail = True
    else:
        _AUX_DIR = os.environ[_AUX_ENV]

        # Check auxil dir is OK...
        for _ in _AUX_FILES:
            if not os.path.exists('%s/%s' % (_AUX_DIR, _)):
                print('Error: %s not in auxil directory!' % _)
                fail = True

    # CALDB check
    if _CALDB_ENV not in os.environ:
        print('Please set the %s environment variable first.' % _CALDB_ENV)
        fail = True
    elif _CALDB_CONF_ENV not in os.environ:
        print('Please set the %s environment variable first.' % _CALDB_CONF_ENV)
        fail = True
    else:
        cal = caldb.CalDB(os.environ[_CALDB_ENV])
        if cal._index is None:
            fail = True

    if fail is True:
        print('Initial check did not pass!')

    return not fail

def block():
    """
    Return True if conf.check() did not pass, and prints an info message.
    """
    if _PASSCHECK is False:
        print('Errors encountered when checking setup for nuskybgd.')
        print('Please fix the problems before continuing.')
        return True
    else:
        return False
