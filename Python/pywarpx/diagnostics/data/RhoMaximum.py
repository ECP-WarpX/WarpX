import sys
import os
import re 
import pandas as pd 
import collections
from warnings import warn

from BaseReader import DataReader

class RhoMaximumData(DataReader):
    """
    Reader for the Rho reduced diagnostic.
    """

    def __init__(self, run_directory, file_prefix='RhoMaximum'):
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
        
        warn("The RhoMaximum diagnostic does not include photon species!!!")

    def get_nonpho_species_names(self):
        '''Returns the names of the non-photon species involved in the diagnostic.'''
        # deduces names of species from the header by considering only the relevant entries 
        # finds the string by locating the units 
        nonpho_species_names = []
        data_file_path = self.get_data_path()       
        header = pd.read_csv(data_file_path, header=0, sep=" ", nrows=0)
        colnames = header.columns.values
        for cn in colnames:
            result=re.search("]max_"+"(.+?)"+"_\|rho\|_lev0", cn)
            if result is not None:
                nonpho_species_names.append(result.group(1))   
        return nonpho_species_names

    def get_nonpho_nspecies(self):
        return len(self.get_nonpho_species_names())


    def get_nlev(self):
        """Returns the number of mesh refinement levels in the simulation."""

        nlev = int((self.get_ncols() - 2.)/(2.+self.get_nonpho_nspecies()))
        return nlev

    def get_valid_args(self):
        species = self.get_nonpho_species_names()
        nlev = self.get_nlev() 
        nspc = self.get_nonpho_nspecies()
        valid_args = [st for n in range(nlev) for st in ['max_rho_lev'+str(n), 'min_rho_lev'+str(n)]]
        valid_args = valid_args + ['max_'+sp+'_|rho|_lev'+str(n) for sp in species for n in range(nlev)]
        return valid_args 

    def get_data(self, *args, steps=None, times=None):
        """
        Arguments:
            string: 
                Can be one or more among the valid arguments. 
        Keyword arguments:     
            steps = list or np.array of integers or None (optional): 
                Timesteps at which the desidered output will be returned.
                If equal to None or not specified then all timesteps are given.  
            times = list or np.array of numbers or None (optional):
                The desidered output will be returned at the closest availble times.
                If equal to None or not specified then all timesteps are given  
        Output: 
            pandas dataframe with columns: steps, times, requested data 
        """
        
        valid_args = self.get_valid_args()
        
        # get data using parent class 
        data = DataReader.get_data(self, valid_args, steps, times, *args)
        
        # restrict to the desidered steps or times 
        restricted_data = self.restrict_data(data, steps, times)
        
        # rename columns         
        restricted_data.columns = ['steps', 'times', *[name for name in args]]

        return restricted_data
