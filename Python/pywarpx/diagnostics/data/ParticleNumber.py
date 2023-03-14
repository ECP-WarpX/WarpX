from BaseReader import DataReader
import numpy as np
import pandas as pd


class ParticleNumberData(DataReader):
    """
    Reader for the ParticleNumber reduced diagnostic.
    """

    def __init__(self, run_directory, file_prefix='ParticleNumber'):
        """
        Parameters
        ----------
        run_directory: string
            Path to the run directory of WarpX.
        file_prefix : string
            Name of the file containing the current reduced diagnostic data.
        """

        super().__init__(run_directory)

        self.data_file_prefix = file_prefix
        self.data_file_suffix = '.txt'

        self.has_species = True

    def get_nspecies(self):
        # remove 4 columns: step, time, total_macroparticles, total_weight
        # divide by 2: we save both the total number of macroparticles
        # and the total number of real particles for every species
        return DataReader.get_nspecies(self, subtract=4., divide=2.)

    def get_species_names(self):
        # remove first 3 entries (step, time, total) then select every 1 entry
        return DataReader.get_species_names(self, string='_macroparticles\(\)', start=3, step=1)

    def get_valid_args(self):
        """
        Returns the valid strings to extract the data.
        """

        species_names = self.get_species_names()
        valid_args = ('total_macroparticles', 'total_weight',
                      *(s+'_macroparticles' for s in species_names), *(s+'_weight' for s in species_names))
        return valid_args

    def get_data(self, *args, **kwargs):
        """
        Arguments:
            string:
                Can be one or more among the valid strings.
        Keyword arguments:
            steps = list or np.array of integers, None or 'all' (optional):
                Timesteps at which the desidered output will be returned.
                If equal to None, 'all' or not specified then all timesteps are given.
            times = list or np.array of numbers, None or 'all' (optional):
                The desidered output will be returned at the closest availble times.
                If equal to None, 'all' or not specified then all timesteps are given
        Output:
            pandas dataframe with columns: steps, times, requested data
        """

        data_file_path = self.get_data_path()
        data = np.loadtxt(data_file_path)

        df = pd.read_csv(data_file_path, sep=" ", header=0)
        cols = df.columns

        species_names = self.get_species_names()

        col = []
        data = []

        valid_args = self.get_valid_args()

        for a in args:
            if a not in valid_args:
                raise ValueError('{} is an invalid argument!\n'
                                  'List of valid arguments: \n'
                                  '{}'.format(a, valid_args))

        if('total_macroparticles' in args):
            col.append('[2]total_macroparticles(')

        if('total_weight' in args):
            col.append([c for c in cols if 'total_weight(' in c][0])

        for s in species_names:
            if (s+'_macroparticles' in args):
                col.append([c for c in cols if s+'_macroparticles(' in c][0])
            if (s+'_weight' in args):
                col.append([c for c in cols if s+'_weight(' in c][0])

        if col == []:
            raise ValueError('Could not find any valid column names!\n'
                             'Is the spelling correct?\n'
                             'Did you select any valid data?')

        col = ['#[0]step()','[1]time(s)'] + col
        data = df[col]

        restricted_data = self.restrict_data(data, **kwargs)

        # rename columns
        restricted_data.columns = ['steps', 'times', *[name for name in args]]

        return restricted_data
