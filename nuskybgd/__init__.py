import importlib

importlib.import_module('%s.util' % __name__)
importlib.import_module('%s.caldb' % __name__)
importlib.import_module('%s.rmf' % __name__)
importlib.import_module('%s.conf' % __name__)

from . import conf

conf._PASSCHECK = conf.check()
