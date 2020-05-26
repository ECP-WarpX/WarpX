# Copyright 2017-2020 David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

from .Bucket import Bucket

diagnostics = Bucket('diagnostics', _diagnostics_dict={})

class Diagnostic(Bucket):
    """
    This is the same as a Bucket, but checks that any attributes are always given the same value.
    """
    def add_new_attr_with_check(self, name, value):
        if name.startswith('_'):
            self._localsetattr(name, value)
        else:
            if name in self.argvattrs:
                assert value == self.argvattrs[name], Exception(f'Diagnostic attributes not consistent for {self.instancename}')
            self.argvattrs[name] = value

    def __setattr__(self, name, value):
        self.add_new_attr_with_check(name, value)

