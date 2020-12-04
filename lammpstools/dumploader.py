# modified (mostly dumbed down for my own use) from the
# Pizza.py toolkit:
# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# and also inspired by code written by Majid Moseyabi (University
# of Bristol)

oneline = "Read LAMMPS dump files and extract thermodynamic data"

import numpy as np
from .baseclass import BaseClass

class DumpLoader(BaseClass):

    """
    Class to load per-atom data from a LAMMPS dump file
    into a list of dictionaries (one dictionary for each
    timestep.

    Attributes
    ----------
    filename : string
        Name of the dump file to load data into.
    data : list of dicts
        data[i] is a dictionary which has
        all data from the ith timestep of the dump file
        stored.
    integerquantities : list of strings
        names of all quantites that should be stored as integers
        instead of floats. Default is ['id','type'].
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
    None


    Private Methods
    ---------------
    __init__(self,filename,min_step=0,max_step=1e10,pstatus)

        initialise attributes, BaseClass, and store data

    __read_dump(self,min_step=0, max_step=1e10,pstatus=False)

        workhorse file which actually stores data into the
        self.data attribute
    
    """

    def __init__(self,filename,integerquantities=['id','type'],
                 min_step=0, max_step=1e10,pstatus=False):

        """
        Initialise class attributes and load data.

        Parameters
        ----------
        filename : string
            Name of the dump file to load data into.
        integerquantities : list of strings
            names of all quantites that should be stored as integers
            instead of floats. Default is ['id','type'].
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

        Returns
        -------
        Does not return anything, but load in data from file.
        """
        
        super(DumpLoader,self).__init__()
        self.filename = filename
        self.integerquantities = integerquantities
        self.min_step = min_step
        self.max_step = max_step

        if (self.max_step<self.min_step) :
            raise Exception('Max_step is smaller than min_step.')

        self.pstatus = pstatus
        
        self.data = self.__read_dump()

        return
    
    def __read_dump(self):
        """
        Read data from dump file self.filename and return a list of
        dictionaries with each list item corresponding to a timestep
        of the dump file.

        Example: If you have a simulation with k atoms, and the (per-atom)
        info is dumped into filename n times total, then this function
        returns a list with n items in it. Each item is a dictionary,
        where the 'key' is the name of the measured quantity (e.g. 'x'
        would be the x-coordinate), and value is typically an array (for
        key 'x', the array would be k atoms long).

        Parameters
        ----------
        None
        
        Returns
        -------
        data : list of dicts
            data[i] is a dictionary which has
            all data from the ith timestep of the dump file
            stored.


        """
        
        data=[] # array which will hold all the snapshots
        minstep_flag = False
        read_flag = False
        min_step = self.min_step
        max_step = self.max_step
        pstatus = self.pstatus
        
        
        cnt = 0


        with open(self.filename,'r') as f:
            snapshot = {} # dictionary to store dump data for current timestep
            snapshot['traj'] = self.filename

            while True:

                
                line=f.readline().strip('\n')
                if not line: break
                items = line.split()
                if items[0] == 'ITEM:':
                    if items[1] == 'TIMESTEP':
                        step = int(f.readline().split(' ')[0])
                        if step > max_step :
                            self._max_step_print(max_step,data,pstatus)
                            return data
                        if ( step >= min_step and minstep_flag == False) :
                            self._min_step_print(min_step,pstatus)
                            minstep_flag = True
                        if (step >= min_step and step <= max_step):
                            read_flag = True
                            cnt += 1
                            self._curr_step_print(cnt,step,pstatus)
                        snapshot['step'] =  step
                    if items[1] == 'NUMBER':
                        N = int(f.readline().split()[0])
                        snapshot['N'] =  N

                    if items[1] == 'BOX':
                        line = f.readline().split()
                        box_x = float(line[1]) - float(line[0])
                        line = f.readline().split()
                        box_y = float(line[1]) - float(line[0])
                        line = f.readline().split()
                        box_z = float(line[1]) - float(line[0])
                        snapshot['z_size'] = box_z
                        snapshot['box'] = np.array([box_x, box_y])

                    if items[1] == 'ATOMS':

                        header = items[2:]
                        x = np.zeros( (N, len(header)) )

                        for i in range(N):
                            line = f.readline().strip('\n').split()
                            for j in range(len(header)):
                                x[ int(line[0])-1, j] = float(line[j])
                                
                        self._set_snap_atom_vals(header,x,snapshot,
                                                 self.integerquantities)

                        if (read_flag): data.append(snapshot.copy())
                        snapshot = {}
                        snapshot['traj'] = self.filename

        return data

    

if __name__ == "__main__":
    
    fname = 'test_data/dump_dummy.lammpstrj'
    
    d = DumpLoader(fname,max_step=50)

    print(d.averages(d.data,sortingcolumn='id'))
