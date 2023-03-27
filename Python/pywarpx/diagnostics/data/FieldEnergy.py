import sys
import os
import re 
import collections

from BaseReader import DataReader

class FieldEnergyData(DataReader):
    """
    Reader for the FieldEnergy reduced diagnostic
    """

    def __init__(self, red_diags_dir="./diags/reducedfiles", separator=" ", fname='FieldEnergy.txt'):
        """
        Parameters
        ----------
        red_diags_dir: string
            Path to the reduced diagnostic directory, default is './diags/reducedfiles'.
        fname: string
            Name of the reduced diagnostic file.
        separator: string
            The separator between row values in the output file, default is a whitespace " ".
        """         
        
        super().__init__(red_diags_dir, separator)

        self.fname = fname
        if self.fname is None:
            raise ValueError('The diagnostic file name cannot be None!')
        
        self.has_species = False 
    
    def get_nlev(self):
        """Returns the number of mesh refinement levels in the simulation."""
 
        return int((self.get_ncols() - 2.)/3.)          
        
    def get_valid_args(self):
        nlev = self.get_nlev() 
        valid_args = [s for n in range(nlev) for s in ['total_lev'+str(n), 'E_lev'+str(n), 'B_lev'+str(n)]]
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
        
        
