"""Array type which can be appended to in an efficient way.
"""

import numpy
import types
# Class which allows an appendable array.
# DPG 8/19/99


class AppendableArray:
    """
    Creates an array which can be appended to in an efficient manner. The object
    keeps an internal array which is bigger than the actual data. When new data is
    appended, it is only copied to fill in the extra space. More space is only
    allocated when that space fills up. This saves both the allocation time, and
    the time to copy the existing data to the new space. This code was written by
    David P. Grote for Warp(8/19/99).

    Arguments:
        initlen (int): The initial size of the array. The most efficiency is gained
            when initlen is close to the total expected length.
        unitshape (tuple): The appendable unit can be an array. This gives the shape of
            the unit. The full shape of the array then would be
            [n]+unitshape, where n is the number of units appended.
        typecode (str): Typecode of the array. Uses the same default as the standard
            array creation routines.
        autobump (int): The size of the increment used when additional extra space
            is needed.
        initunit (np.ndarray): When given, the unitshape and the typecode are taken from
            it. Also, this unit is make the first unit in the array.
        aggressivebumping (float): Whenever new space is added, autobump will be
            increased by this factor in an attempt to
            minimize the number of reallocations. Its new
            value will be the max of the size of the appended data
            and this times the old autobump size. A good value is
            1.5 - this can greatly reduce the amount of
            rallocation without a significant amount of wasted
            space.

    Create an instance like so
    >>> a = AppendableArray(initlen=100,typecode='d')
    Append a single unit like this
    >>> a.append(7.)
    or multiple units like this
    >>> a.append(numpy.ones(5,'d'))
    The data can be obtained by directly indexing the instance, like this
    >>> print a[:4]
    [ 7., 1., 1., 1.,]
    will give the first four number appended
    Other methods include len, data, setautobump, cleardata, reshape
    """
    def __init__(self,initlen=1,unitshape=None,typecode=None,autobump=100,
                 initunit=None,aggressivebumping=1.5):
        if typecode is None: typecode = numpy.zeros(1).dtype.char
        self._maxlen = initlen
        if initunit is None:
            self._typecode = typecode
            self._unitshape = unitshape
        else:
            # --- Get typecode and unitshape from initunit
            if isinstance(initunit, numpy.ndarray):
                self._typecode = initunit.dtype.char
                self._unitshape = initunit.shape
            elif isinstance(initunit, int):
                self._typecode = 'i'
            else:
                self._typecode = 'd'
                self._unitshape = None

        self._datalen = 0
        self._autobump = autobump
        self.aggressivebumping = aggressivebumping
        self._allocatearray()
        if initunit is not None:
            self.append(initunit)

    def checkautobumpsize(self,deltalen):
        # --- A factor of 1.5 gives nearly the same amount of savings as 2,
        # --- but doesn't waste quite as much space.
        newautobump = max(deltalen,int(self.aggressivebumping*self.getautobump()))
        self.setautobump(newautobump)

    def _extend(self,deltalen):
        # --- Only increase of the size of the array if the extra space fills up
        if len(self) + deltalen > self._maxlen:
            self.checkautobumpsize(deltalen)
            self._maxlen = self._maxlen + max(deltalen,self._autobump)
            a = self._array[:len(self),...] + 0
            self._allocatearray()
            self._array[:len(self),...] = a

    def _allocatearray(self):
        if self._unitshape is None:
            self._array = numpy.zeros(self._maxlen,self._typecode)
        else:
            self._array = numpy.zeros([self._maxlen]+list(self._unitshape),self._typecode)

    def append(self,data):
        if self._unitshape is None:
            # --- If data is just a scalar, then set length to one. Otherwise
            # --- get length of data to add.
            try:
                lendata = len(data)
            except (TypeError,IndexError):
                lendata = 1
        else:
            # --- Data must be an array in this case.
            # --- If the shape of data is the same as the original shape,
            # --- then only one unit is added. Otherwise, get the number
            # --- of units to add. The length is always added to the first
            # --- dimension.
            if len(data.shape) == len(self._unitshape): lendata = 1
            else:                                       lendata = data.shape[0]
        self._extend(lendata)
        newlen = self._datalen + lendata
        self._array[self._datalen:newlen,...] = data
        self._datalen = newlen

    def data(self):
        """
        Return the data.
        """
        return self._array[:len(self),...]

    def setautobump(self,a):
        """
        Set the autobump attribute to the value specified.
        """
        self._autobump = a

    def getautobump(self):
        """
        Get the autobump attribute
        """
        return self._autobump

    def cleardata(self):
        """
        Reset the array so it has a length of zero.
        """
        self._datalen = 0

    def resetdata(self, data):
        """
        Resets the data to be the input values - all of the original data is thrown
        away. The unit shape of data must be the same.
        """
        self.cleardata()
        self.append(data)

    def compressdata(self, data, whereidx):
        """
        Takes in data and copies into preallocated array using boolean indexing

        Arguments:
            data (array):  raw data to be copied
            whereidx (boolean array): must be same size as first index in data
        """

        # Make sure the input data is an array and the dimensions match the unitshape
        if not isinstance(data, numpy.ndarray):
            raise Exception("data: must be a numpy array")

        # Make sure where is a boolean array
        if ((not isinstance(whereidx, numpy.ndarray))
            or (whereidx.dtype != bool)
            or (len(whereidx) != data.shape[0])):
            raise Exception("whereidx must be a numpy array of booleans, with"\
                            + " the same size as the first dimension of data")

        self.cleardata()

        # Make sure there is enough room in the array
        datalen = whereidx.sum()
        self._extend(datalen)
        self._datalen = datalen

        # Copy compressed data into array
        outview = self._array[:self._datalen,...]
        numpy.compress(whereidx, data, axis=0, out=outview)

    def takedata(self, data, idx):
        """
        Takes in data and copies into preallocated array using fancy indexing

        Arguments:
            data (array):  raw data to be copied
            idx (int array): must be equal or less in size as first index in data
        """

        # Make sure the input data is an array and the dimensions match the unitshape
        if not isinstance(data, numpy.ndarray):
            raise Exception("data: must be a numpy array")

        # Make sure where is an int array
        if ((not isinstance(idx, numpy.ndarray))
            or (idx.dtype != int)):
            raise Exception("idx must be a numpy array of ints")

        self.cleardata()

        # Make sure there is enough room in the array
        datalen = len(idx)
        self._extend(datalen)
        self._datalen = datalen

        # Copy compressed data into array
        outview = self._array[:self._datalen,...]
        numpy.take(data, idx, axis=0, out=outview)


    def unitshape(self):
        if self._unitshape is None:
            return (1,)
        else:
            return self._unitshape

    def reshape(self,newunitshape):
        """
        Change the shape of the appendable unit. Can only be used if a unitshape was
        specified on creation.

        Arguments
            newunitshape (tuple): must have the same number of dimensions as the original
                unitshape
        """
        assert self._unitshape is not None,\
               'Only an array with a specified unitshape can be reshaped'
        assert len(newunitshape) == len(self._unitshape),\
               ('New unitshape must have the same number of dimensions as original ' +
                'unitshape')
        # --- Save old data
        oldunitshape = self._unitshape
        oldarray = self._array
        # --- Create new array
        self._unitshape = newunitshape
        self._allocatearray()
        # --- Copy data from old to new
        ii = [None] + list(numpy.minimum(oldunitshape,newunitshape))
        ss = tuple(map(slice,ii))
        self._array[ss] = oldarray[ss]

    def __len__(self):
        return self._datalen

    def __getitem__(self,key):
        return self.data()[key]

    def __setitem__(self,key,value):
        self.data()[key] = value
