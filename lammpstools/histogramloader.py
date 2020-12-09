# modified (mostly dumbed down for my own use) from the
# Pizza.py toolkit:
# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# and also inspired by code written by Majid Moseyabi (University
# of Bristol)

oneline = "Read LAMMPS dump files and extract thermodynamic data"

import numpy as np
from .baseclass import BaseClass

class HistogramLoader(BaseClass):

    """
    Class to load non-local data (e.g. binned data such
    as the radial distribution function) into a list of
    dictionaries (one dictionary for each timestep.

    Attributes
    ----------
    filename : string
        Name of the dump histogram to load data from.
    data : list of dicts
        data[i] is a dictionary which has
        all data from the ith timestep of the histogram file
        stored.
    integerquantities : list of strings
        names of all quantites that should be stored as integers
        instead of floats. Default is ['Row'].
    min_step : int (optional)
        Any timesteps in the dump file which are less than min_step are
        ignored (useful for initial transients). Default is 0.
    max_step : int
        Any timesteps in the dump file which are greater than max_step are
        ignored (useful if dump file is massive and you don't want to
        store it all). Default is 1e10.
    pstatus : bool (optional)
        Print where the function is at in terms of how much data it has
        stored (kind of like a verbose option).

    Public Methods
    --------------
    averages(self)

        Compute averages over timesteps per row. Derived from
        averages method in base class (see docs there).

    Constructor
    -----------
    __init__(self,filename,min_step=0,max_step=1e10,pstatus)

        initialise attributes and store data

    Private Methods
    ---------------
    __read_histograms(self,min_step=0, max_step=1e10,pstatus=False)

        workhorse file which actually stores data into the
        self.data attribute
    

    """

    def __init__(self,fname,skiplines=2,min_step=0,
                 max_step=1e10,pstatus=False,
                 integerquantities=['Row']):
        """
        Initialise class attributes and load data.

        Parameters
        ----------
        fname : string
            Name of the dump histogram to load data from.
        skiplines : int (optional)
            Number of lines to skip in the file before the header file (e.g.
            if third line of file reads '# Row c_rdf[1] c_rdf[2] c_rdf[3]',
            then you would want skiplines = 2
        min_step : int (optional)
            Any timesteps in the dump file which are less than min_step are
            ignored (useful for initial transients). Default is 0.
        max_step : int
            Any timesteps in the dump file which are greater than max_step are
            ignored (useful if dump file is massive and you don't want to
            store it all). Default is 1e10.
        pstatus : bool (optional)
            Print where the function is at in terms of how much data it has
            stored (kind of like a verbose option).
        integerquantities : list of strings
            names of all quantites that should be stored as integers
            instead of floats. Default is ['Row'].

        Returns
        -------
        Does not return anything, but loads in data.
        """
        
        self.fname = fname
        self.min_step = min_step
        self.max_step = max_step
        
        self.integerquantities = integerquantities


        if (self.max_step<self.min_step) :
            raise Exception('Max_step is smaller than min_step.')

        self.pstatus = pstatus
        
        self.data = self.__read_histograms(skiplines)

        return

    def averages(self):
        """
        Compute averages over timesteps per row. Derived from
        averages method in base class (see docs there).
        """

        return super(HistogramLoader,self).averages(self.data)
    
    def __read_histograms(self,skiplines):
        """
        Read data from histogram file self.fname and return a list of
        dictionaries with each list item corresponding to a timestep
        of the histogram file.

        Example: If you have a simulation with k atoms, and the (per-atom)
        info is dumped into fname n times total, then this function
        returns a list with n items in it. Each item is a dictionary,
        where the 'key' is the name of the measured quantity (e.g.
        'c_rdf[1]' would be the radial distance in the radial
        distribution function), and value is typically an array (for
        key 'c_rdf[1]', the array would be nbins long).

        Parameters
        ----------
        skiplines : int
            Number of lines to skip in the file before the header file (e.g.
            if third line of file reads '# Row c_rdf[1] c_rdf[2] c_rdf[3]',
            then you would want skiplines = 2
            
        
        Returns
        -------
        data : list of dicts
            data[i] is a dictionary which has
            all data from the ith timestep of the histogram file
            stored.


        """
        
        data=[] # array which will hold all the snapshots

        minstep_flag = False
        read_flag = False
        min_step = self.min_step
        max_step = self.max_step
        pstatus = self.pstatus
        
        
        cnt = 0


        with open(self.fname,'r') as f:
            
            count = 0
            
            while count < skiplines:
                line=f.readline()
                count += 1

            items = f.readline().strip('\n').split()
            header = items[1:]
                
            while True:
                
                snapshot = {} # store binned data at current timestep
                snapshot['traj'] = self.fname
                
                line = f.readline().strip('\n')
                if not line:
                    break
                
                items = line.split()
                
                timestep = int(items[0])
                if timestep < min_step:
                    continue
                elif timestep > max_step:
                    self._max_step_print(max_step,data,pstatus)
                    break
                else:
                    if not minstep_flag:
                        self._min_step_print(min_step,pstatus)
                    minstep_flag=True
                    cnt += 1
                    self._curr_step_print(cnt,timestep,pstatus)

                
                snapshot['timestep'] = timestep
                snapshot['nbins'] = int(items[1])

                N = snapshot['nbins']
                n = len(header)
                x = np.zeros( (N, n) )

                for i in range(N):
                    
                    line = f.readline().strip('\n').split()
                    
                    for j in range(n):
                        
                        x[ int(line[0])-1, j] = float(line[j])
                self._set_snap_atom_vals(header,x,snapshot,
                                         self.integerquantities)
                data.append(snapshot.copy())
        return data
            

if __name__ == "__main__":
    
    fname = 'test_data/histograms_vsfp_0.001_0.0_1.0_3.0_1.rdf'
    
    rdf = HistogramLoader(fname,pstatus=True)

    avdata, stddata = rdf.averages()

    print(avdata)
