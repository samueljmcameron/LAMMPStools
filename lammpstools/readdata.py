import numpy as np


class ReadData():

    """

    Read in a LAMMPS formatted data file. See
    https://lammps.sandia.gov/doc/read_data.html for
    documentation on what a LAMMPS data file consists of.

    The point of this is to be used as a way to test new
    LAMMPS code. E.g. to test a new compute, compute
    the result with the lammps code and then also with
    python, both cases the data file (the latter needing
    this class).

    Attributes
    ----------

    headers : dict
        dictionary with keys being any header keywords
        in the file (e.g. 'atom', 'improper types', etc)
        and values being the corresponding value.
    atoms : dict
        dictionary with keys being atom-property keywords
        listed in the file (e.g. 'Atoms', 'Ellipsoids', etc)
        and values being themselves dictionaries holding
        relevant properties (e.g. atom_props['Atoms'] would
        be a dictionary containing atom-ID, atom-type, x y z
        coords, and more depending on the atom_style).

    Public Methods
    --------------
    
    get_atom_style(self)
        Get the atom style of the class.

    """

    def __init__(self, fname, atom_style, skiprows=1):

        self._atom_style = atom_style
        self.headers,self.atoms = self._read_data(fname,skiprows)

        return

    def get_atom_style(self):

        """
        Get the atom style of the class.

        Returns
        -------
        atom_style : string
            atom style as set by the user when the class was
            instantiated.
         
        """

        return self._atom_style
    
    def _read_data(self,fname,skiprows):

        """
        Workhorse method which actually reads data from input file
        and organises it.

        Currently only supports "Atoms" keywords for the body of the
        file.

        Parameters
        ----------
        fname : string
            file name of the data file to be read
        skiprows : int
            number of rows to skip in data file before
            starting to read it

        Returns
        -------
        header_dict : dict
            dictionary with keys being any header keywords
            in the file (e.g. 'atom', 'improper types', etc)
            and values being the corresponding value.
        atom : dict
            dictionary with keys being atom-property keywords
            listed in the file (e.g. 'Atoms', 'Ellipsoids', etc)
            and values being themselves dictionaries holding
            relevant properties (e.g. atom_props['Atoms'] would
            be a dictionary containing atom-ID, atom-type, x y z
            coords, and more depending on the atom_style).

        """

        
        with open(fname,"r") as datfile:

            for i in range(skiprows):
                datfile.readline()

            header_dict = {}
            line = self._next_nonempty_line(datfile)
            headerentry,headerout = self._check_if_header(line)
            while headerentry != None:

                header_dict[headerentry] = headerout
                line = self._next_nonempty_line(datfile)
                headerentry,headerout = self._check_if_header(line)


            atom = {}
            Natoms = header_dict["atoms"]
            while (line != ''):

                linelist,comments = self._line2list(line)

                if linelist[0] == "Atoms":

                    self._Atoms_warning(comments)
                    
                    atom["Atoms"] = self._mkAtoms(datfile,Natoms)

                line = self._next_nonempty_line(datfile)

        return header_dict,atom

    def _Atoms_warning(self,comments):

        """
        Check whether read_data file "Atoms" keyword is consistent
        with the user's prescribed atom_style when instantiating. E.g.
        if the data file says 'Atoms # sphere' but atom_style variable
        is set to 'atomic', then this function will print a warning.

        Parameters
        ----------
        comments : list of strings
            comment words after 'Atoms' keyword, if it exists. E.g.
            ['sphere', 'hello'] in 'Atoms # sphere hello'. Could
            be empty.

        Returns
        -------
        Does not return anything.
        """
        

        atom_styles = ["angle","atomic","body","bond","charge","dipole",
                       "dpd","edpd","electron","ellipsoid","full","line",
                       "mdpd","mesont","molecular","peri","smd","sph",
                       "sphere","spin","tdpd","template","tri",
                       "wavepacket"]


        if (len(comments) > 0 and comments[0] in atom_styles and 
            comments[0] != self._atom_style):
            
            print("Warning! You have set atom_style = '{self._atom_style}', "
                  "but read_data indicates that the 'Atoms' keyword is "
                  "followed by '# {comments[0]}' indicating conflicting "
                  "styles.")

        return


    def _mkAtoms(self,datfile,Natoms):

        """
        Initialise and fill dictionary which holds per-atom
        properties (from keyword 'Atoms' in data file).

        Currently only supports atom_style 'atomic'.

        Parameters
        ----------
        datfile : read-able file (from e.g. open(fname,'r'))
            LAMMPS data file, with current position at the 
            'Atoms' line.
        Natoms : int
            number of atoms
        
        Returns
        -------
        Atoms : dict
            dictionary containing per-atom properties as specified
            by the 'Atoms' body keyword.

        """
        
        Atoms = {}


        if self._atom_style == 'atomic':

            Atoms['atom-ID'] = np.empty([Natoms],float)
            Atoms['atom-type'] = np.empty([Natoms],float)
            Atoms['x'] = np.empty([Natoms],float)
            Atoms['y'] = np.empty([Natoms],float)
            Atoms['z'] = np.empty([Natoms],float)

            
            for i in range(Natoms):
                line = self._next_nonempty_line(datfile)
                linelist,comments = self._line2list(line)
                Atoms['atom-ID'][i] = int(linelist[0])
                Atoms['atom-type'][i] = int(linelist[1])
                Atoms['x'][i] = float(linelist[2])
                Atoms['y'][i] = float(linelist[3])
                Atoms['z'][i] = float(linelist[4])
            

        else:
            raise ValueError("atom_style must be 'atomic' (for now).")

        return Atoms

    
    def _line2list(self,line):
        """
        Convert a raw line (e.g. a string as read by readline())
        into two lists of words, one containing pre-comment
        words and another containing post-comment words (the comment
        delimiter '#' is not included in either list).

        Parameters
        ----------
        line : string
            line to be converted to a list

        Returns
        -------
        pre : list of strings
            list of words preceding '#' (can be empty).
        post : list of strings
            list of words preceded by '#' (can be empty).

        """
        # check if comment is present in line

        if '#' in line:

            cindex = line.index('#')
            # add space between '#' and next char if necessary

            if line[cindex+1] != ' ':
                line = line[:cindex+1] + ' ' + line[cindex+1:]


        rawout = line.strip('\n').split(' ')
        out = []
        for word in rawout:
            if word != '':
                out.append(word)

        if '#' in out:
            ic = out.index('#')
        else:
            ic = len(out)

        return out[:ic],out[ic+1:]
    
    def _next_nonempty_line(self,datfile):

        """
        Find next non-empty line in a file.

        Parameters
        ----------
        datfile : read-able file (from e.g. open(fname,'r'))
            LAMMPS data file.

        Returns
        -------
        line : string
            Next non-empty line in the file, or empty string if there
            are no more non-empty lines.

        """

        line = datfile.readline()
        while (line == '\n'):
            line = datfile.readline()

        return line


    def _check_wordnum_header(self,wordlist,entry,nwords):

        """
        Given a header entry, check if the number of words is what
        it should be. If it's not, raise a value error.

        Parameters
        ----------
        wordlist : list of strings
            words in the header entry line
        entry : string
            name of header entry line (e.g. 'atoms',
            'improper types', etc)
        nwords : int
            expected number of words

        Returns
        -------
        Does not return anything.

        """


        if len(wordlist) != nwords:
            raise ValueError(f"Too many words in '{entry}' header line.")

        return



    def _check_if_header(self,line):

        """
        Given a line, process it if its a LAMMPS data file header entry.

        Parameters
        ----------
        line : string
            string which may or may not be a LAMMPS data file header entry.

        Returns
        -------
        entry : string or None
            Name of header entry keyword (if line is truly a header line),
            otherwise it returns None
        value : int or string or list or tuple (whatever really) or None
            Value of header entry (preceding keyword) if line is truly a
            header line, otherwise returns None.



        """


        header_entries = ["atoms", "bonds", "angles", "dihedrals", "impropers",
                          "atom types", "bond types", "angle types", "dihedral types",
                          "improper types", "ellipsoids", "lines", "triangles",
                          "bodies", "xlo xhi", "ylo yhi", "zlo zhi", "xy xz yz"]


        linelist,comments = self._line2list(line)




        if linelist[-1] == 'types':
            entry = linelist[-2] + ' ' + linelist[-1]

            if entry not in header_entries:
                return None,None

            # test for errors
            self._check_wordnum_header(linelist,entry,3)

            # get value preceding header entry
            out = int(linelist[0])

        elif linelist[-1][-2:] == 'hi':
            entry = linelist[-2] + ' ' + linelist[-1]


            if entry not in header_entries:
                return None,None


            # test for errors
            self._check_wordnum_header(linelist,entry,4)

            # get value preceding header entry
            lo = float(linelist[0])
            hi = float(linelist[1])

            out = (lo,hi)

        elif linelist[-1] == 'yz':
            entry = linelist[-3] + ' ' + linelist[-2] + ' ' + linelist[-1]

            if entry not in header_entries:
                return None,None

            # test for errors
            self._check_wordnum_header(linelist,entry,6)

            # get value preceding header entry
            xy = float(linelist[0])
            xz = float(linelist[1])
            yz = float(linelist[2])


            out = (xy,xz,yz)

        else:
            entry = linelist[-1]

            if entry not in header_entries:
                return None,None

            # test for errors
            self._check_wordnum_header(linelist,entry,2)

            out = int(linelist[0])


        return entry,out

        
        

    
