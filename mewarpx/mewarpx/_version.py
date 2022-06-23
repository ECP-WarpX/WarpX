# One and only one place to store the version info
# https://stackoverflow.com/questions/458550/standard-way-to-embed-version-into-python-package
__version_info__ = (8, 4, 1)
__version__ = '.'.join([str(x) for x in __version_info__])

# One and only one place to store the Physics version
__physics_version__ = 2
