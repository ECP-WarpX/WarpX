import os
import re 
import numpy as np
import pandas as pd 
import collections

from BaseReader import DataReader

class ParticleMomentumData(DataReader):
    """
    Reader for the ParticleMomentum reduced diagnostic.
    """

    def __init__(self, run_directory, file_prefix='ParticleMomentum'):
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
        # remove 8 columns: step, time, total_x, total_y, total_z, total_mean_x, total_mean_y, total_mean_z
        # divide by 6: for each species we save both the total and mean momentum along the 3 coordinates 
        return DataReader.get_nspecies(self, subtract=8., divide=2.*3.)
                     
    def get_species_names(self):
        return DataReader.get_species_names(self, string='_x\(kg\*m\/s\)', start=5, step=3)
    
    def get_valid_args(self):
        """
        Returns the valid strings to extract the data. 
        """    
        
        species_names = self.get_species_names()
        valid_args = ('total_x', 'total_y', 'total_z',  
                      'total_mean_x', 'total_mean_y', 'total_mean_z',
                      *(s+'_x' for s in species_names), *(s+'_y' for s in species_names), *(s+'_z' for s in species_names), 
                      *(s+'_mean_x' for s in species_names), *(s+'_mean_y' for s in species_names), *(s+'_mean_z' for s in species_names))        
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
               
        if('total_x' in args):
            col.append('[2]total_x(kg*m/s)')
        if('total_y' in args):
            col.append('[3]total_y(kg*m/s)')
        if('total_z' in args):
            col.append('[4]total_z(kg*m/s)')

        if('total_mean_x' in args):
            col.append([c for c in cols if 'total_mean_x(kg*m/s)' in c][0])
        if('total_mean_y' in args):
            col.append([c for c in cols if 'total_mean_y(kg*m/s)' in c][0])            
        if('total_mean_z' in args):
            col.append([c for c in cols if 'total_mean_z(kg*m/s)' in c][0])
        
        for s in species_names:
            if (s+'_x' in args):
                col.append([c for c in cols if s+'_x(kg*m/s)' in c][0])  
            if (s+'_y' in args):
                col.append([c for c in cols if s+'_y(kg*m/s)' in c][0])  
            if (s+'_z' in args):
                col.append([c for c in cols if s+'_z(kg*m/s)' in c][0])                  
            if (s+'_mean_x' in args):
                col.append([c for c in cols if s+'_mean_x(kg*m/s)' in c][0])   
            if (s+'_mean_y' in args):
                col.append([c for c in cols if s+'_mean_y(kg*m/s)' in c][0])                         
            if (s+'_mean_z' in args):
                col.append([c for c in cols if s+'_mean_z(kg*m/s)' in c][0])  
                                
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
