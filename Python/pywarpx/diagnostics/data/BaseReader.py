import collections
import os
import re

import numpy as np
import pandas as pd


class DataReader():
    """
    Base class that all data readers should inherit from.
    Works only in Cartesian geometry.
    """

    def __init__(self, run_directory, file_prefix=None):
        """
        Parameters
        ----------
        run_directory: string
            Path to the run directory of WarpX.
        """
        if run_directory is None:
            raise ValueError('The run_directory parameter can not be None!')

        self.run_directory = run_directory
        self.data_file_prefix = None
        self.data_file_suffix = None

        self.has_species = None

    def get_data_path(self, subdir='reducedfiles'):
        """
        Returns the path to the data file.
        Parameters
        ----------
        subdir: string
            Name of reducedfiles directory, default is 'reducedfiles'.
        Returns
        -------
        A string with a path to the current diagnostic (assumed unique).
        """

        sim_output_dir = os.path.join(self.run_directory, 'diags/'+subdir+'/')
        if not os.path.isdir(sim_output_dir):
            raise IOError('The simulation directory does not exist inside '
                          'path:\n  {}\n'
                          'Did you set the proper path to the run directory?\n'
                          'Did the simulation already run?'
                          .format(self.run_directory))

        data_file_path = os.path.join(
            sim_output_dir,
            self.data_file_prefix +
            self.data_file_suffix
        )
        if not os.path.isfile(data_file_path):
            raise IOError('The file {} does not exist.\n'
                          'Did the simulation already run?'
                          .format(data_file_path))

        return data_file_path

    def get_ncols(self):
        ''' Returns the number of columns of the reduced diagnostic.'''
        data_file_path = self.get_data_path()
        header=pd.read_csv(data_file_path, header=0, sep=" ", nrows=0)
        ncols = len(header.columns)
        return ncols

    def get_nspecies(self, subtract=None, divide=None):
        '''Returns the number of species involved in the diagnostic.'''
        # deduces number of species from number of columns
        if self.has_species == True:
            ncols = self.get_ncols()
            # remove <substract> columns and divide to avoid counting multiple entries
            nspecies = int((ncols - subtract)/divide)
            return nspecies
        else:
            return 0

    def get_species_names(self, string=None, start=None, step=None):
        '''Returns the names of the species involved in the diagnostic.'''
        # deduces names of species from the header by considering only the relevant entries
        # finds the string by locating the units
        if self.has_species == True:
            species_names = []
            data_file_path = self.get_data_path()
            nspecies = self.get_nspecies()

            header = pd.read_csv(data_file_path, header=0, sep=" ", nrows=0)
            colnames = header.columns.values
            for i in range(nspecies):
                result=re.search('](.+?)'+string,colnames[start+i*step])
                species_names.append(result.group(1))
            return species_names
        else:
            raise ValueError('The current diagnostics does not involve any species!')

    def load_data(self):
        data_file_path = self.get_data_path()
        df = pd.read_csv(data_file_path, sep=" ", header=0)
        return df

    def get_steps(self):
        data_file_path = self.get_data_path()
        df = pd.read_csv(data_file_path, sep=" ", usecols=['#[0]step()'])
        steps = np.unique(df.values)
        return steps

    def get_times(self):
        data_file_path = self.get_data_path()
        df = pd.read_csv(data_file_path, sep=" ", usecols=['[1]time(s)'])
        times = np.unique(df.values)
        return times

    def get_data(self, valid_args, steps=None, times=None, *args):

        df = self.load_data()
        all_cols = df.columns

        col = []
        data = []

        for a in args:
            if a not in valid_args:
                raise ValueError('{} is an invalid argument!\n'
                                  'list of valid arguments: \n'
                                  '{}'.format(a, valid_args))
            else:
                col.append([c for c in all_cols if a in c][0])

        if col == []:
            raise ValueError('Could not find any valid column names!\n'
                             'Is the spelling correct?\n'
                             'Did you select any valid data?')

        col = ['#[0]step()','[1]time(s)'] + col
        data = df[col]

        return data


    def restrict_data(self, data, steps=None, times=None): # check if works for all diags (e.g. probes, field reduction)
        """
        Restricts the data to the requested steps or times.
        Arguments:
            pandas dataframe:
                Extracted data.
        Keyword arguments:
            steps = list or np.array of integers or None (optional):
                Timesteps at which the desidered output will be returned.
                If equal to None or not specified then all timesteps are given.
            times = list or np.array of numbers or None (optional):
                The desidered output will be returned at the closest available time.
                If equal to None or not specified then all timesteps are given.
        Output:
            pandas dataframe with columns: steps, times, extracted data
        """

        all_times = self.get_times()
        all_steps = self.get_steps()

        if (steps is not None) and (times is not None):
            raise ValueError('Please provide only one between steps and times!')
        elif steps is not None:
            # make steps iterable if they are not
            if not isinstance(steps, collections.abc.Iterable):
                steps = np.array([steps])
        elif times is not None:
            if not isinstance(times, collections.abc.Iterable):
                times = np.array([times])
            steps = []
            for t in times:
                i = np.argmin((t-all_times)**2)
                steps = np.append(steps, all_steps[i])
        else:
            steps = all_steps

        # verify that requested iterations exist
        if not set(steps).issubset(all_steps):
            raise IndexError('Selected step {} is not available!\n'
                             'List of available steps: \n'
                             '{}'.format(steps, all_steps))

        data = data.loc[data['#[0]step()'].isin(steps)]
        # reset index
        data = data.reset_index(drop=True)

        return data
