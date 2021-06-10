# One and only one place to store the version info
# https://stackoverflow.com/questions/458550/standard-way-to-embed-version-into-python-package
__version_info__ = (0, 1, 0)
__version__ = '.'.join([str(x) for x in __version_info__])
