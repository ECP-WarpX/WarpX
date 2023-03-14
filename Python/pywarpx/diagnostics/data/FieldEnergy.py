import sys
import os
import re 
import numpy as np
import pandas as pd 
import collections

from BaseReader import DataReader

class FieldEnergyData(DataReader):
    """
    Reader for the FieldEnergy reduced diagnostic
    """

    def __init__(self, run_directory, file_prefix='FieldEnergy'):
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
        
        self.has_species = False 
    
    def get_nlev(self):
        """Returns the number of mesh refinement levels in the simulation."""
 
        return int((self.get_ncols() - 2.)/3.)          
        
    def get_valid_args(self):
        nlev = self.get_nlev() 
        valid_args = [s for n in range(nlev) for s in ['total_lev'+str(n), 'E_lev'+str(n), 'B_lev'+str(n)]]
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
        
        col = []
        data = []
        
        nlev = self.get_nlev()      
        valid_args = self.get_valid_args()

        for a in args:
            if a not in valid_args:
                raise ValueError('{} is an invalid argument!\n' 
                                  'List of valid arguments: \n'
                                  '{}'.format(a, valid_args))
        for n in range(nlev):
            if('total_lev'+str(n) in args):
                col.append([c for c in cols if 'total_lev'+str(n)+'(J' in c][0])
            if('E_lev'+str(n) in args):
                col.append([c for c in cols if 'E_lev'+str(n)+'(J' in c][0])
            if('B_lev'+str(n) in args):
                col.append([c for c in cols if 'B_lev'+str(n)+'(J' in c][0])
                
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
