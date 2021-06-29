"""
Utility functions for mewarpx.
"""
import inspect
import os

from pywarpx import geometry

# http://stackoverflow.com/questions/50499/in-python-how-do-i-get-the-path-and-name-of-the-file-t
mewarpx_dir = os.path.dirname(os.path.abspath(
    inspect.getfile(inspect.currentframe())))


def init_libwarpx(ndim, rz):
    """_libwarpx requires the geometry be set before importing.

    This complicates a lot of our code if we need to delay importing it - so
    instead we import it here.

    Very Bad Things will happen if ndim and rz here are different than is
    used in the rest of the simulation!

    Arguments:
        ndim (int): Number of dimensions. Ignored for RZ.
        rz (bool): True for RZ simulations, else False.
    """
    geometry.coord_sys = 1 if rz else 0
    geometry.prob_lo = [0]*ndim
    import pywarpx._libwarpx
    # This just quiets linters like pyflakes by using the otherwise-unused
    # variable
    assert pywarpx._libwarpx
