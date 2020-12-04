
import numpy as np

class BaseClass(object):

    """
    Base class used for non-thermodynamic data output
    (e.g. dump files or radial distribution function
    files).


    Public Methods
    --------------
    averages(self,data,sortingcolumn='Row')

        Compute average values and standard deviations over all time
        steps of data by row.


    Constructor
    ---------------
    __init__(self)

    Public Methods
    --------------

    Protected Methods
    ---------------

    _max_step_print(self,max_step,data,pstatus)
        
        Print a message to stdout that the final time-step in the
        file has been reached.

    _min_step_print(self,min_step,pstatus)

        Print a message to stdout that the first time-step in the
        file has been reached.

    _curr_step_print(self,cnt,step,pstatus)
        Print a message to stdout that the current time-step in the
        file is being read.

    _set_snap_atom_vals(self,header,x,snapshot)

         Insert values of the file into snapshot dictionary.


    """

    def __init__(self):

        return
    

    def _max_step_print(self,max_step,data,pstatus):
        """
        Print a message to stdout that the final time-step in the
        file has been reached.

        Enabled if pstatus is true.
        """
        
        if pstatus:
            print(f'max_step reached {max_step}')
            print(f"Last TIMESTEP {data[-1]['step']}")
        return

    def _min_step_print(self,min_step,pstatus):
        """
        Print a message to stdout that the first time-step in the
        file has been reached.

        Enabled if pstatus is true.
        """
        
        if pstatus:
            print(f'First TIMESTEP reached {min_step}')
        return

    def _curr_step_print(self,cnt,step,pstatus):
        """
        Print a message to stdout that the current time-step in the
        file is being read.

        Enabled if pstatus is true.
        """
        
        if pstatus:
            print(f'{cnt}: reading TIMESTEP {step}')
        return

    def _set_snap_atom_vals(self,header,x,snapshot,
                             integerquantities):
        """
        Insert values of the file into snapshot dictionary.

        Parameters
        ----------
        header : list
            list of strings containing names of all the quantities
            measured in the x (next parameter), e.g. 'x', 'xu', etc.

        x : np.ndarray
            contains all quantities measured in the timestep.
            E.g. if the the file only outputs 'x' and 'y'
            positions of 20 atoms (or 20 bins, etc), then x
            will be an array of shape 2 by 20.

        snapshot : dict
            dictionary where data of current timestep gets stored,
            e.g. snapshot[header[i]] =  x[:,i]

        """
        

        for j, key in enumerate(header):
            if key in integerquantities : 
                snapshot[key] = np.array( [int(p) for p in x[:,j]]  ) 
            else:
                snapshot[key] = x[:,j]

        return
    
    def averages(self,data,sortingcolumn='Row'):
        """
        Compute average values and standard deviations over all time
        steps of data by row.

        Parameters
        ----------
        data : list of dicts
            data[i] is a dictionary which has
            all data from the ith timestep of a simulation.
            Note that all dicts in data must have the same
            keys.
        sortingcolumn : string (optional)
            Name of key in the list of dicts which should
            be the same for each item in the list. Required
            to confirm that it makes sense to take the average
            (see notes below). For example, in LAMMPS dump file,
            the 'id' key will indicate this (and will only 
            satisfy the above condition if the dump file is
            sorted by atom id). Default value is 'Row',
            because this method is most useful in histogram
            type data (e.g. compute rdf output from LAMMPS).

        Returns
        -------
        avdata : dict
            Average values of each key in the list of dicts
            from data.
        stddata : dict
            Standard deviation (as computed by numpy)
            of the values of each key in the list of dicts
            from data.

        Notes
        -----
        The main use of this function is on histogram like
        files, where all the time steps have binned values
        at the same e.g. radial distances.

        In contrast, when used on a LAMMPS dump file, this
        function will compute the average values of
        quantities for each individual atom over all time
        steps (e.g. if atom 1 has x=4 at t=0, x=8 at t=1,
        and atom 2 has x =-2 at t = 0 and x=-3 at t=1, then
        this function would return {'x' : np.array([6,-2.5])}.
        Therefore, it only makes sense to call this function
        on dump files if they are sorted. 

        """
        

        # dictionary of average values
        avdata = {}
        stddata = {}

        d0 = data[0]
        Ntotal = len(data)
        nrows = len(d0[sortingcolumn])
            
        mega_dict = {}

        errmsg = ("Cannot average unsorted data by row. Check that "
                  "all timesteps have same column array in key "
                  f"{sortingcolumn}.")
        
        for i,snap in enumerate(data):

            for key,value in snap.items():
                if (isinstance(value, np.ndarray)
                    and len(value)==nrows):
                    if i == 0:
                        if key in self.integerquantities:
                            mega_dict[key] = np.empty([nrows,Ntotal],
                                                      int)
                        else:
                            mega_dict[key] = np.empty([nrows,Ntotal],
                                                      float)
                    mega_dict[key][:,i] = value
                        
                    if key == sortingcolumn:
                        if i > 0:
                            if not np.isclose(mega_dict[key][:,i],
                                              mega_dict[key][:,i-1]).all():
                                
                                raise RuntimeError(errmsg)

        for key, value in mega_dict.items():

            avdata[key] = np.mean(value,axis=1)
            stddata[key] = np.std(value,axis=1)

        return avdata,stddata
