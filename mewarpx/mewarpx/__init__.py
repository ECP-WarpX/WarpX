"""mewarpx package"""
import logging
import sys

from mewarpx.utils_store import mwxlogging

from ._version import __version__, __version_info__

__author__ = 'Modern Electron <peter.scherpelz@modernelectron.com>'
__all__ = []

# This runs logging.basicConfig() if it has not been run already. This may not
# be the most flexible idea, but if we aren't using logging otherwise, it
# should make things easier. See
# https://lists.gt.net/python/python/863016 and
# https://docs.python.org/2/howto/logging.html and
# https://stackoverflow.com/questions/44188270/no-handlers-could-be-found-for-logger
logger = logging.getLogger(__name__)
logger.debug("Initializing logging for mewarpx")

h = mwxlogging.MEHandler(sys.stdout)
f = mwxlogging.MEFilter()

h.addFilter(f)
logger.addHandler(h)

logger.setLevel(logging.INFO)
