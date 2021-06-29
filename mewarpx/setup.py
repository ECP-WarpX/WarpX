""" Runs pip to install this and requirements. See Makefile for commands."""
# Template for this package from
# https://github.com/kragniz/cookiecutter-pypackage-minimal
import os
import setuptools

requires = [
    "numpy",
    "scipy",
]

extras = {
    "tests": ["pytest", "pytest-cov", "pytest-xdist"],
    # "docs": ["sphinx", "mock", "sphinx_rtd_theme"],
    # "aws": ["awscli", "boto3", "s3fs"]
}
# http://stackoverflow.com/questions/19096155/setuptools-and-pip-choice-of-minimal-and-complete-install
complete_extras = [i for sublist in extras.values() for i in sublist]
extras["complete"] = complete_extras

# Load version file - see
# https://stackoverflow.com/questions/458550/standard-way-to-embed-version-into-python-package
here = os.path.abspath(os.path.dirname(__file__))
version_file = 'mewarpx/_version.py'
version_file = os.path.join(here, version_file)
if not os.path.isfile(version_file):
    raise IOError(
        "Version file _version.py not found; are you installing from the root "
        "directory?"
    )
with open(version_file, 'r') as vfile:
    exec(vfile.read())

# Remove 'tests' and subdirs as an installed package
setuptools_packages = setuptools.find_packages(
    exclude=["tests", "tests.*", "tests_parallel", "tests_parallel.*"])

setuptools.setup(
    name="mewarpx",
    version=__version__,
    url="https://github.com/ModernElectron/mewarpx",

    author="Modern Electron",
    author_email="peter.scherpelz@modernelectron.com",

    description=(
        "A set of tools for running and postprocessing WarpX simulations with "
        "additional thermionic-related capabilities."
    ),
    long_description=open('README.rst').read(),

    packages=setuptools_packages,

    install_requires=requires,

    # These get installed by installing extras
    extras_require=extras,

    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    # entry_points={
    #     'console_scripts': [
    #     ]
    # },
    #
    # test_suite='tests',
    # tests_require=[],
)
