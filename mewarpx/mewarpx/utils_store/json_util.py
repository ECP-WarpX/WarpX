from collections import OrderedDict
import datetime
import json

import numpy as np
import periodictable

from mewarpx import diags, emission


class MEWarpXEncoder(json.JSONEncoder):
    """A JSON encoder for MEWarpX objects."""

    def default(self, object):
        """Default operations if JSON library fails to serialize an object."""
        if isinstance(object, datetime.datetime):
            return object.strftime((r"%m/%d/%Y_%H:%M:%S:%f"))

        elif isinstance(object, OrderedDict):
            return object

        elif isinstance(object, (bytes, bytearray)):
            return object.decode("utf-8")

        # only include name of element
        elif isinstance(object, (periodictable.core.Isotope,
                                 periodictable.core.Element)):
            return object.name

        # leave out diags to avoid circular serialization error
        elif isinstance(object, diags.WarpXDiagnostic):
            return str(object)

        # remove very long and unneeded attributes
        elif isinstance(object, emission.Emitter):
            obj_dict = object.__dict__.copy()
            obj_dict.pop("xvec", None)
            obj_dict.pop("yvec", None)
            obj_dict.pop("zvec", None)
            obj_dict.pop("contours", None)
            obj_dict.pop("centers", None)
            obj_dict.pop("dvec", None)
            obj_dict.pop("distances", None)
            obj_dict.pop("CDF", None)
            obj_dict.pop("normal", None)
            return obj_dict

        elif isinstance(object, np.ndarray):
            return object.tolist()

        # for generic objects, use object dictionary
        elif hasattr(object, "__dict__"):
            return object.__dict__

        else:
            return super().default(self, object)
