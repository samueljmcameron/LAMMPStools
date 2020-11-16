# modified (mostly dumbed down for my own use) from the
# Pizza.py toolkit:
# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html


oneline = "Read LAMMPS log files and extract thermodynamic data"

import numpy as np

class LogLoader(object):

    """
    Class to load thermodynamic data from a LAMMPS log file into
    either a single dictionary (merging all chunks of data) or a
    list of dictionaries.

    Attributes
    ----------
    data : dict or list of dicts
        Depending on initialisation, i.e. whether user decides
        to merge all thermodynamic chunks together (see __init__#
        docstring for details), either a single dict or a list of
        dicts. The dictionary(s) containing keys which are the
        names of the thermodynamic quantities being measured (e.g.
        'Step', 'Temp', etc.), and the values are the arrays
        data corresponding to each key. So, if the chunks are
        not merged, then data[0]['Temp'] would be an array of all
        the temperatures measured in the first data chunk. If the
        chunks are merged, then data['Temp'] would be an array of
        all the temperatures measured in the log file.

    Public Methods
    --------------
    remove_chunks(self,i)

        Remove one or more of the thermodynamic data chunks from the
        data list. Useful if you want to merge all the data later but
        want to exclude certain chunks.

    merge_chunks(self,remove_duplicates=True)
       
        Merge all of the chunks from the data file into a single
        dictionary.


    Private Methods
    ---------------
    __init__(self,fname,chunk_start = 'Step',chunk_end = 'Loop')

        Initialise all attributes, and apply any chunk removals or
        merging of the chunks.

    __load_data(self)

        Load in all thermodynamic data chunks from log file.             

    __load_lines(self)

        Load in all lines of the log file.
        
    __chunk_lines(self,lines)

        Build an array mask which is true for all data lines
        (i.e. lines with thermo output) and false for all
        other lines.

    __build_header(self,items)

        Construct a new dictionary with keys having the
        values of item in the items list.
    """

    def __init__(self,fname,chunk_start = 'Step',chunk_end = 'Loop'):
        """
        Initialise all parameters, and apply any chunk removals or merging
        of the chunks.

        Parameters
        ----------
        fname : string
            The name of the log file.
        chunk_start : string (optional)
            A string which signifies the start of a thermodynamic data
            chunk within the log file. Almost always this will be 'Step'.
        chunk_end : string (optional)
            A string which signifies the end of a thermodynamic data
            chunk within the log file. Almost always this will be 'Loop'.
        """

        self._fname = fname
        self._chunk_start = chunk_start
        self._chunk_end = chunk_end
        
        self.data = self.__load_data()

        return
    # --------------------------------------------------------------------

    def remove_chunks(self,i):
        """
        Remove one of the thermodynamic data chunks from the the data
        list. Useful if you want to merge all the data later but
        want to exclude a specific chunk.

        Parameters
        ----------
        i : float or array
            Index(es) of data list component to be removed.

        Returns
        -------
        Does not return anything, but overwrites the existing
        attribute data (list of dicts), excluding the chunks
        just removed.
        """

        data = []

        if type(i) is np.array or type(i) is list:
            count = 0
            i = sorted(i)
            
            for j,d in enumerate(self.data):
                if j != i[count]:
                    data.append(d)
                else:
                    count += 1
                    
        elif type(i) is int:
                        
            for j,d in enumerate(self.data):
                if j != i:
                    data.append(d)

        else:
            raise Exception('parameter i must be either an int, '
                            'list, or np.array')

        self.data = data
        
        return


    # --------------------------------------------------------------------

    def merge_chunks(self,remove_duplicates=True):
        """
        Merge all of the chunks from the data file into a single
        dictionary.

        Parameters
        ----------
        remove_duplicates : boolean (optional)
            If True, remove any duplicate data points (defined as any
            two data points which have the same 'Step' value), keeping
            the data point which comes later in the file.

        Returns
        -------
        dataset : dict
            Dictionary of thermodynamic data.
        """
        
        dataset = dict(self.data[0])
        
        for j,ddict in enumerate(self.data[1:]):
            
            for key,value in ddict.items():
                try:
                    val = dataset[key]
                except KeyError:

                    raise Exception('Unable to merge thermo chunks: '
                                    f'{key} variable in thermo '
                                    'chunk {j} is not measured in '
                                    'the previous thermo chunk.')
                dataset[key] = np.concatenate((val,
                                               ddict[key]))

        if remove_duplicates:
            steps = np.copy(dataset['Step'])
            
            dmask = steps[1:]-steps[:-1]
            dmask = np.array(np.concatenate(([True],dmask)),bool)

            for key in dataset:
                dataset[key] = dataset[key][dmask]
            
            assert (dataset['Step'] == np.unique(steps)).all()

            
        return dataset
    
    # --------------------------------------------------------------------
    
    def __load_data(self):
        """
        Load in all thermodynamic data chunks from log file.

        Returns
        -------
        data : list
            List where each component is a dictionary containing
            thermodynamic data (the keys of each dictionary correspond
            to the names of the different thermodynamic observables,
            the values are arrays with the actual data).
             
        """

        data = []

        lines = self.__load_lines()
        
        master_masks,header_list = self.__chunk_lines(lines)

        for j,chunk in enumerate(header_list):
            
            d_lines = np.array(lines)[master_masks[j]]
            
            dt = np.array([d.split() for d in d_lines])

            for key,value in chunk.items():

                if key == 'Step':
                    chunk[key] = np.array(dt[:,value],int)
                else:
                    chunk[key] = np.array(dt[:,value],float)
            data.append(chunk)
                                    
        return data

    
    # --------------------------------------------------------------------
    
    def __load_lines(self):
        """
        Load in all lines of the log file.
        
        Returns
        -------
        lines : list
            List of all lines in the log file. 
        """
        with open(self._fname) as f:
            
            lines = [line.strip('\n') for line in f.readlines()
                     if line.strip()]

        return lines

    
    # --------------------------------------------------------------------

    def __chunk_lines(self,lines):
        """
        Build an array mask which is true for all data lines
        (i.e. lines with thermo output) and false for all
        other lines.

        """
        
        chunk_start = self._chunk_start
        chunk_end = self._chunk_end
        header_list = []
        master_masks = []

        cnt = 0
        while cnt < len(lines):
            
            line = lines[cnt]
            items = line.split()
            
            # if chunk of thermo outputs starts here
            
            if items[0] == chunk_start:

                # construct header containing all thermo variables
                # that will be output, with keys being variable name
                # and value being the column of the variable
                headerdict = self.__build_header(items)

                # construct mask which will select only lines of
                # current thermo output chunk
                mask = np.zeros_like(lines,bool)

                # move to next line, since current line has variable
                # names
                cnt += 1

                line = lines[cnt]
                items = line.split()
                
                # while still reading numbers

                while items[0] != chunk_end:
                    # if reading a slurm file, sometimes MPI
                    # produces messages in the middle of data
                    # chunks about processors, so check to
                    # ensure the line in the data chunk is
                    # actually a data line, which is why the
                    # following try is necessary


                    try:
                        int(items[0])
                    except:
                        print("weird interruption in data stating: "
                              f"{' '.join(items)}")
                    else:
                        mask[cnt] = True

                    cnt += 1
                    line = lines[cnt]
                    items = line.split()

                # add header and masks to lists
                header_list.append(headerdict)
                master_masks.append(mask)

            else:
                cnt+= 1
                
        return master_masks,header_list

    def __build_header(self,items):
        """
        Construct a new dictionary with keys having the
        values of item in the items list.

        Parameters
        ----------
        items : list
            List of objects (usually strings) to be used as
            keynames in the dictionary.

        Returns
        -------
        newdict : dict
            A dictionary with keys from the items list.
        """
        
        newdict = {}

        for j,item in enumerate(items):
            newdict[item] = j

        return newdict
            

if __name__=='__main__':


    prefix = './'
    
    fname = 'log_100_0.5_3.0_1.lammps.log'

    ll = LogLoader(fname)

    print(ll.remove_chunks.__doc__)
    print(ll.merge_chunks.__doc__)
