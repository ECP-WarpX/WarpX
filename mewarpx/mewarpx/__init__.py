"""mewarpx package"""
import logging

from ._version import __version__, __version_info__

__author__ = 'Modern Electron <peter.scherpelz@modernelectron.com>'
__all__ = []

# This runs logging.basicConfig() if it has not been run already. This may not
# be the most flexible idea, but if we aren't using logging otherwise, it
# should make things easier. See
# https://lists.gt.net/python/python/863016 and
# https://docs.python.org/2/howto/logging.html and
# https://stackoverflow.com/questions/44188270/no-handlers-could-be-found-for-logger
logging.debug("Initializing logging for mewarpx")
