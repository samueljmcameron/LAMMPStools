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
    __init__(self,filename,integerquantities=['id','type'],
             min_step=0,max_step=1e10,pstatus=False)

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

    def clustering(self,snapshot,clusterlabel):
        """ Store sorted cluster information into snapshot dictionary.

        Creates five new keys in snapshot dict:

        'cluster_sizes' value is a dict with keys labeling the
        cluster number and values counting number of atoms in the
        cluster.

        'N_clusters' value is an integer (number of unique clusters).

        'lammps2sorted_id' value is a dictionary, keys are labels of
        original cluster label (as outputted by dump), and values are
        labels of sorted ordering of cluster sizes (i.e. going from
        0 to N_clusters -1).

        'sorted2lammps_id' value is similar to 'lammps2sorted_id'
        above, except keys and values are switched.

        'new_c_c1' is the labels of cluster sizes, maintaining
        the same order as 'c_c1', except that the cluster labels
        have been changed such that the labels of atoms from the
        largest cluster in 'c_c1' (which are completely random)
        are represented as 0s in 'new_c_c1', and the labels of
        atoms from the smallest cluster in 'c_c1' will be labelled
        as N_clusters-1.
        e.g. if snapshot['c_c1'] = [1,2,3,4,4,4,5,1,6], then
        snapshot['new_c_c1'] = [1,2,3,0,0,0,5,6]. Note that
        ordering is random for those clusters with the same
        number of atoms, so the above is equivalent to
        snapshot['new_c_c1'] = [1,6,5,0,0,0,2,3], since there
        are four atoms which are of cluster size one.

        """
        
        cluster_sizes = self.__count_items(snapshot[clusterlabel])
        snapshot['cluster_sizes'] = cluster_sizes

        snapshot['N_clusters'] = len(cluster_sizes)
        
        # compute list of cluster LABELS, going from largest 
        # cluster to smallest cluster. Note this is NOT a list
        # of cluster sizes.
        sorted_key = sorted(cluster_sizes, key=cluster_sizes.get,
                            reverse=True)

        lammps2sorted_id = { key:c_id  for c_id, key
                             in enumerate(sorted_key)   }
        sorted2lammps_id = { c_id:key  for c_id, key
                             in enumerate(sorted_key)   }
        snapshot['lammps2sorted_id'] = lammps2sorted_id
        snapshot['sorted2lammps_id'] = sorted2lammps_id
        snapshot['new_c_c1'] = [lammps2sorted_id[lmp_id]
                                for lmp_id in snapshot[clusterlabel]]

        return

    def __count_items(self,lst):
        """ Count number of times a label occurs in a list.

        Input lst is the list with labels. 

        Returns a dictionary, with keys being all unique items from
        lst and values being the count of how many times those items
        occur in  lst.

        An example form of this list from a dump file would be the
        cluster (c_c1) row. If the dump entry had 4 atoms, and the
        first and third were in the same cluster, lst = [1,2,1,3]
        and this function would return {'1' : 2, '2' : 1, '3' : 1}.

        """
        
        my_count = {}
        for val in lst:
            my_count[val] = (my_count.get(val, 0)) + 1
        return my_count

    

if __name__ == "__main__":
    
    fname = 'test_data/dump_dummy.lammpstrj'
    
    d = DumpLoader(fname,max_step=50)

    print(d.averages(d.data,sortingcolumn='id'))
