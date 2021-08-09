from builtins import object
import collections
import math

import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage


def concat_crop_timeseries(timeseries_list,
                           step_begin=None,
                           step_end=None,
                           keys=None,
                           dt=None,
                           debug=False):
    """Combine or crop timeseries object(s).

    Note:
        dt must be a multiple of all individual timeseries' dt. Timeseries that
        do not share the final dt will be resampled with default smoothing.

    Arguments:
        timeseries_list (list of Timeseries): List of objects to concatenate.
            If empty, return None.
        step_begin (int): Step to begin the array with, cropping timeseries if
            needed. If None (default), take all timesteps from inputs.
        step_end (int): Step to end the array with + 1, cropping timeseries if
            needed. If None (default), take all timesteps from inputs.
        dt (float): The timestep to use. If None, use the maximum dt among the
            series. In all cases, this must be a multiple of each timeseries
            dt.
        keys (list of str): If specified, concatenate only these keys and ignore
            any others held by the objects.
        debug (bool): If True, print out lots of diagnostic info. Default
            False. These print statements are commented out by default to
            prevent frequent if statement checks; so uncomment to enable debug
            functionality.

    Returns:
        timeseries (Timeseries): The concatenated, cropped, Timeseries object.
    """
    if len(timeseries_list) < 1:
        # debug_print(debug, "No timeseries to concatenate")
        return None

    # Get dt. If it doesn't work for resampling, error will be thrown during
    # resampling.
    if dt is None:
        dt = max([x.dt for x in timeseries_list])
        # debug_print(debug, "dt auto-calculated as ", dt)

    ts_resampled = []
    for ts in timeseries_list:
        # debug_print(
        #     debug,
        #     "ts has {} elements".format(ts.step_end - ts.step_begin)
        # )
        ts_resampled.append(ts.resample(new_dt=dt, smooth=True))
        # debug_print(
        #     debug,
        #     "resampled ts has step_begin {}, step_end {}".format(
        #         ts_resampled[-1].step_begin, ts_resampled[-1].step_end)
        # )

    # Get step begin/end
    if step_begin is None:
        step_begin = min(x.step_begin for x in ts_resampled)
    # debug_print(debug, "step_begin is ", step_begin)

    if step_end is None:
        step_end = max(x.step_end for x in ts_resampled)
    # debug_print(debug, "step_end is ", step_end)

    # Create new object
    final_ts = Timeseries(step_begin=step_begin, step_end=step_end,
                          dt=dt)

    # Iterate through Timeseries and each key.
    if keys is None:
        keys = set([y for x in timeseries_list
                    for y in list(x._array_dict.keys())])
        # debug_print(debug, "keys", keys)

    for key in keys:
        final_ts._array_dict[key] = np.zeros(step_end - step_begin)

    for ts in ts_resampled:

        # Figure out which parts of arrays to transfer by finding begin and
        # ending indices for each.
        if ts.step_begin >= step_begin:
            ts_beginidx = 0
            self_beginidx = ts.step_begin - step_begin
        else:
            ts_beginidx = step_begin - ts.step_begin
            self_beginidx = 0

        if ts.step_end >= step_end:
            ts_endidx = step_end - ts.step_begin
            self_endidx = step_end - step_begin
        else:
            ts_endidx = ts.step_end - ts.step_begin
            self_endidx = ts.step_end - step_begin

        # debug_print(debug, "Found ts idxs", ts_beginidx, ts_endidx)
        # debug_print(debug, "Found self idxs", self_beginidx, self_endidx)

        # In cases where a timeseries falls entirely outside the new range, we
        # skip it completely. This avoids a bug where you can get eg
        # ts_beginidx = 0 and ts_endidx = -1; the concatenation will fail.
        if ts_endidx <= ts_beginidx:
            continue

        for key in keys:
            # debug_print(debug, "Concat", key)
            if key not in ts._array_dict:
                continue

            # debug_print(debug, "Data",
            #             ts._array_dict[key][ts_beginidx:ts_endidx])
            final_ts._array_dict[key][self_beginidx:self_endidx] += (
                ts._array_dict[key][ts_beginidx:ts_endidx]
            )

    return final_ts


class Timeseries(object):

    """Hold a regularly-spaced array of values that we can then manipulate."""

    def __init__(self, step_begin, step_end, dt, array_dict=None):
        """Initialize timeseries object with the arrays for it to hold.

        Arguments:
            step_begin (int): First timestep stored in the arrays
            step_end (int): Last timestep stored in the arrays
            dt (float): Time increment between steps in seconds
            array_dict (dict of np.ndarray): A dictionary with keys that are
                strings, and values that are one-dimensional arrays with length
                (step_end - step_begin)/stride, containing the value of that
                string's quantity at each step = step_begin + i*stride. If None,
                initialize an empty dictionary that can be used later. All
                quantities should already be rates for resampling to work (eg
                current, not charge, at a given timestep).
        """
        self.step_begin = int(round(step_begin))
        self.step_end = int(round(step_end))
        self.dt = dt
        self._array_dict = array_dict
        if self._array_dict is None:
            self._array_dict = collections.OrderedDict()

        self._check_input()

    def _check_input(self):
        """Perform error checking on data."""
        if self.step_end <= self.step_begin:
            raise ValueError("Timeseries must have step_end > step_begin")

        if self.dt <= 0.:
            raise ValueError("dt must be greater than 0.")

        for key, val in self._array_dict.items():
            if len(val) != self.n_elements:
                raise ValueError(
                    f"All arrays must have correct number of elements {self.n_elements} "
                    f"based on step values, but array {key} has {len(val)} elements."
                )

    @property
    def n_elements(self):
        """A property so that it automatically adapts to resampling."""
        return self.step_end - self.step_begin

    def keys(self):
        """A user accessible list of keys."""
        return list(self._array_dict.keys())

    def set_array(self, key, array):
        """Add or set an array in the dictionary.

        Arguments:
            key (str): Key to add
            array (np.ndarray): Array to add
        """
        self._array_dict[key] = array
        self._check_input()

    def get_timeseries_by_key(self, key, include_times=True, default=None):
        """Get timeseries, including an array of times if requested.

        Arguments:
            key (str): Key to look up.
            include_times (bool): If True, return an n_elements x 2 array that
                has time in seconds in first column; values in second. If False,
                return an n_elements array with values only.
            default (None or float): Value to fill timeseries array with if the
                given key does not exist in the dictionary or if the dictionary
                is empty. If None, AttributeError will be raised for an invalid
                query.
        """
        if self._array_dict is None:
            array_vals = default
        else:
            array_vals = self._array_dict.get(key, default)

        if array_vals is None:
            raise AttributeError(
                '%s not found in the dictionary or the dictionary is None.'%key
            )

        if not include_times:
            return array_vals

        timeseries_array = np.empty((self.n_elements, 2))
        timeseries_array[:, 0] = (
            np.arange(self.step_begin, self.step_end) * self.dt
        )
        timeseries_array[:, 1] = array_vals

        return timeseries_array

    def get_averagevalue_by_key(self, key, default=None):
        """Get the rate value of a timeseries.

        Arguments:
            key (str): Key to look up.
            default (float): Value to return if the given key does not exist in
                the dictionary or if the dictionary is empty. If None,
                AttributeError will be raised for an invalid query.

        Returns:
            value (float): Average value
        """
        if self._array_dict is None:
            array_vals = default
        else:
            array_vals = self._array_dict.get(key, default)

        if array_vals is None:
            raise AttributeError(
                '%s not found in the dictionary or the dictionary is None.'%key
            )

        return np.mean(array_vals)

    def resample(self, new_dt, inplace=False, smooth=True):
        """Resample the time series at a longer timescale.

        Arguments:
            new_dt (float): Must be a multiple of existing dt!
            inplace (bool): If True, modify this Timeseries. If False, return a
                new one. Note that if resample is the identity operation, the
                current Timeseries is returned regardless of this setting.
                Default False.
            smooth (bool): If True, first apply Gaussian smoothing with sigma =
                0.5*new_dt to the timeseries before sampling from the
                timeseries. Default True. Note smoothing is never applied if
                new_dt is equal to current dt. If resampling by 500x or more,
                smooth using mean values instead.

        Returns:
            timeseries (Timeseries): Object with new dt.
        """
        if new_dt < self.dt:
            raise ValueError(
                "new_dt must be greater or equal to existing dt."
            )

        dt_factor = int(round(new_dt / self.dt))

        if not np.isclose(new_dt / self.dt, dt_factor):
            raise ValueError(
                "new_dt must be a multiple of existing dt."
            )

        if dt_factor == 1:
            return self

        step_begin = self.step_begin // dt_factor
        step_begin_remainder = self.step_begin % dt_factor

        if step_begin_remainder > 0:
            step_begin += 1
            idx_offset = dt_factor - step_begin_remainder
        else:
            idx_offset = 0

        # This calculation makes sure that when there are not an evenly
        # divisible number of elements, we get the correct step_end - otherwise
        # there are off-by-one errors with step_end one too small. Eg
        # step_begin=1, step_end=11, dt_factor=2 should give new_step_begin=1,
        # new_step_end=6, and sample elements 2, 4, 6, 8, 10 from original.
        n_elements = int(math.ceil(
            (self.n_elements - idx_offset)
            / float(dt_factor)
        ))
        step_end = step_begin + n_elements

        if inplace:
            target = self
            self.dt = new_dt
            self.step_begin = step_begin
            self.step_end = step_end
        else:
            target = Timeseries(step_begin=step_begin,
                                step_end=step_end,
                                dt=new_dt,
                                )

        for key, val in self._array_dict.items():
            if smooth:
                if dt_factor < 500:
                    target._array_dict[key] = (
                        self.gauss_smooth(
                            val,
                            dt_factor/2.
                        )[idx_offset::dt_factor]
                    )
                else:
                    target._array_dict[key] = (
                        self.mean_smooth(val, dt_factor, idx_offset)
                    )
            else:
                target._array_dict[key] = val[idx_offset::dt_factor]

        target._check_input()

        return target

    @staticmethod
    def gauss_smooth(array, sigma):
        """Apply Gaussian smoothing to the given array.
        This function is based on minerva.util.gaussian_smoothing, but since it
        operates in a given way on 1D arrays only, with no resampling, it's
        re-implemented in simplified form here.

        Arguments:
            array (np.ndarray): The 1D array to smooth.
            sigma (float): The sigma, in units of steps, to use for smoothing.

        Returns:
            smoothed_array (np.ndarray): 1D array with smoothing applied.
        """
        # Both truncate=4 and mode=reflect are the default values for this
        # function.
        smoothed_array = ndimage.gaussian_filter(
            array,
            sigma=sigma,
            truncate=4.0,
            mode='reflect'
        )

        return smoothed_array

    @staticmethod
    def mean_smooth(array, dt_factor, idx_offset=0):
        """Apply smoothing and resampling to the given array.
        This function is a simplified smoothing scheme appropriate for cases
        where many data points will be combined into a single point. Smoothing
        is done by simply averaging the points to be combined.

        Arguments:
            array (np.ndarray): The 1D array to smooth.
            dt_factor (int): The units of steps to use for smoothing. If 1 or
                smaller, then no smoothing is performed and the original array
                is returned as is. If len(array) or larger, then the entire
                array is smoothed down to a single point.
            idx_offset (int): Smoothing calculation is performed on the slice
                array[idx_offset:]. Entries 0:idx_offset are ignored in the
                input and not returned in the output. Defaults to 0.

        Returns:
            smoothed_array (np.ndarray): 1D array with smoothing applied.
        """
        array_length = len(array) - idx_offset

        # Handle edge cases
        if dt_factor < 1:
            dt_factor = 1
        if dt_factor > array_length:
            dt_factor = array_length

        new_array_length = int(np.ceil(array_length / dt_factor))
        full_intervals = int(array_length / dt_factor)
        end_offset = idx_offset + array_length - array_length % dt_factor

        smoothed_array = np.zeros(new_array_length)

        smoothed_array[:full_intervals] = np.mean(
            np.split(array[idx_offset:end_offset], full_intervals),
            axis=1
        )

        if new_array_length > full_intervals:
            smoothed_array[-1] = np.mean(
                array[idx_offset + full_intervals * dt_factor:]
            )

        return smoothed_array


class TimeseriesPlot(object):

    """Handle plotting of arbitrary timeseries."""

    # Start with default params
    default_timeseries_params = {
        "xlabel": r'Time ($\mu$s)',
        "ylabel": 'Timeseries output',
        "title": 'Timeseries plot',
        "labelsize": 24,
        "legendsize": 14,
        "titlesize": 32,
        "alpha": 1.0,
        "yfactor": 1.0,
        "legend_loc": 'best',
        "legend_bbox_to_anchor": None,
    }

    def __init__(self, array_list, ax=None, **kwargs):
        r"""Plot fluxes throughout the simulation.

        Arguments:
            array_list (list of tuples of (name, timeseries_array)):
                timeseries_arrays are those returned by
                ``Timeseries.get_timeseries_by_key(key, include_times=True)``.
                name is the label for the plot legend for that timeseries. List
                order determines plotting order. Smoothing should be performed
                before being passed to this function.
            ax (matplotlib.Axes): If specified, plot on this axes object. If
                unspecified, get current axes.
            xlabel (string): abscissa label. Must include 'ns' or '$\mu$s' to
                specify units in this string.
            ylabel (string): ordinate label
            title (string)
            labelsize (int): default 24
            legendsize (int): default 14
            titlesize (int): default 32
            alpha (float): Transparency of lines to avoid obscuring lines, if
                desired.  Default 1 (disabled).
            yfactor (float): Multiply y-axis by this value. Eg for correcting
                the effective area in currents.
            legend_loc (str): 'loc' argument in legend call
            legend_bbox_to_anchor (tuple of float): bbox_to_anchor argument in
                legend call
        """
        self.timeseries_params = self.default_timeseries_params.copy()
        self.timeseries_params.update(kwargs)

        if ax is None:
            ax = plt.gca()

        # determine scaling factor for time-axis
        if 'ns' in self.timeseries_params["xlabel"]:
            scaling_factor = 1e9
        elif r'$\mu$s' in self.timeseries_params["xlabel"]:
            scaling_factor = 1e6
        else:
            raise ValueError(
                r"Unrecognized time in xlabel. Specify 'ns' or '$\mu$s'."
            )

        for (name, ts_array) in array_list:
            ax.plot(scaling_factor*ts_array[:, 0],
                    self.timeseries_params["yfactor"]*ts_array[:, 1],
                    alpha=self.timeseries_params["alpha"], label=name)

        if kwargs.get('legend', True):
            ax.legend(fontsize=self.timeseries_params["legendsize"],
                      loc=self.timeseries_params["legend_loc"],
                      bbox_to_anchor=self.timeseries_params[
                          "legend_bbox_to_anchor"])
        ax.set_xlabel(self.timeseries_params["xlabel"],
                      fontsize=self.timeseries_params["labelsize"])
        ax.set_ylabel(self.timeseries_params["ylabel"],
                      fontsize=self.timeseries_params["labelsize"])
        ax.set_title(self.timeseries_params["title"],
                     fontsize=self.timeseries_params["titlesize"])
