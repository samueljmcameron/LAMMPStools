import numpy as np


def PvdToDump(inputpvdname,outputdumpname,boxbounds,vtkscalars=['id','type'],
              vtkvectors={'x':('x','y','z'),'F':('Fx','Fy','Fz'),'ux':('ux','uy','uz')},
              integerquantities=['id','type'],dt=None,headerorder = None,
              vtpprefix=None):


    """
    Convert a pvd file (which references a series of vtp files) into something like
    a lammps dump file. The assumed format of the pvd file is that it contains DataSets
    with either:
        a) .vtp files with filename style f'somefilename_{timestep}.vtp'
        b) .pvtp files with filename style f'somefilename_{timestep}.pvtp'
    Furthermore, the vtp or pvtp files must be in the same directory as the pvd file.
    If the pvd file references .pvtp files, then see vtpprefix argument below.
    
    
    Parameters
    ----------
    inputpvdname : string
        Name of the vtk pvd file to read from.
    outputdumpname : string
        Name of the dump file to be written.
    boxbounds : numpy array of shape (3,2)
        Contains the lower and upper limits of the simulation box in x, y, and z directions.
    vtkscalars : list of strings (optional)
        Names of scalar data that should be found in the vtk file. Default is ['id','type'].
    vtkvectors : dictionary (optional)
        Keys are the names of vector data in the vtk file, values are the individual
        component names to be used for the output. Default is 
        {'x':('x','y','z'),'F':('Fx','Fy','Fz'),'ux':('ux','uy','uz')}.
    integerquantities : list of strings (optional)
        Any vtkscalars or values in the vtkvectors dict that are to be treated as integers.
        Default is ['id','type']
    dt : float or None
        If not None, use pvd file 'timestep' information to determine timesteps. If None,
        assume that the series of vtp files ends with the form f"_{timestep}.vtp".
    headerorder : list of strings
        Strings representing the order in which the per atom data will be written.
    vtpprefix : string or None
        If not None, then use this string to find '.vtp' files 
        instead. This assumes .vtp files of the form f'{vtpprefix}_{timestep}.vtp' (without
        any directory prefix needed as the .vtp files should be in the same folder as the
        pvd file).

    Returns
    -------
    Nothing, just writes a dump file.
    """

    # determine if pvd and vtp files are in a different directory
    folder = inputpvdname[:inputpvdname.rfind("/")+1]
    
    with open(inputpvdname,"r") as pvdfile:
        with open(outputdumpname,"w") as dumpfile:    
            for line in pvdfile:
                if line[:8] != "<DataSet":
                    continue
                substr = line[line.find("timestep=\"")+10:]
                time = float(substr[:substr.find("\"")])

                substr = line[line.find("file=\"")+6:]
                vtkfname = folder + substr[:substr.find("\"")]
                
                if dt != None:
                    tstep = round(time/dt)
                else:
                    tstep = int(vtkfname[vtkfname.rfind("_")+1:vtkfname.find(".")])

                if vtpprefix != None:
                    vtkfname = folder + f"{vtpprefix}_{tstep}.vtp"

                vtkToDump(vtkfname,dumpfile,tstep,boxbounds,vtkscalars=vtkscalars,
                          vtkvectors=vtkvectors,integerquantities=integerquantities,
                          headerorder=headerorder)

    return
        


def vtkToDump(vtkfname,dumpfile,timestep,boxbounds,vtkscalars=['id','type'],
              vtkvectors={'x':('x','y','z'),'F':('Fx','Fy','Fz'),'ux':('ux','uy','uz')},
              integerquantities=['id','type'],headerorder = None):
    """
    Convert a vtk file with per atom data (a .vtp file) to something like a lammps dump file.

    Parameters
    ----------

    vtkfname : string
        Name of the vtk file to read from.
    dumpfile : writable file object
        File to write the data to.
    timestep : int
        Timestep of vtk file.
    boxbounds : numpy array of shape (3,2)
        Contains the lower and upper limits of the simulation box in x, y, and z directions.
    vtkscalars : list of strings (optional)
        Names of scalar data that should be found in the vtk file. Default is ['id','type'].
    vtkvectors : dictionary (optional)
        Keys are the names of vector data in the vtk file, values are the individual
        component names to be used for the output. Default is 
        {'x':('x','y','z'),'F':('Fx','Fy','Fz'),'ux':('ux','uy','uz')}.
    integerquantities : list of strings (optional)
        Any vtkscalars or values in the vtkvectors dict that are to be treated as integers.
        Default is ['id','type']
    headerorder : list of strings
        Strings representing the order in which the per atom data will be written.
    
    """

    print(f"Writing timestep {timestep} to dump file.")
    data,natoms = vtkToDict(vtkfname,vtkscalars=vtkscalars,
                            vtkvectors=vtkvectors,
                            integerquantities=integerquantities)

                
    dumpfile.write(f"ITEM: TIMESTEP\n{timestep}\n")
    dumpfile.write(f"ITEM: NUMBER OF ATOMS\n{natoms}\n")
    dumpfile.write(f"ITEM: BOX BOUNDS pp pp pp\n")
    for i in range(3):
        dumpfile.write(f"{boxbounds[i,0]} {boxbounds[i,1]}\n")

    dumpfile.write("ITEM: ATOMS ")
    if headerorder == None:
        headerorder = data.keys()

    for header in data.keys():
        if header not in headerorder:
            raise ValueError(f"headerorder specified is missing {header}")
        
    for header in headerorder:
        dumpfile.write(header + " ")
    dumpfile.write("\n")
    for atom in range(natoms):
        for header in headerorder:
            dumpfile.write(str(data[header][atom]) + " ")
        dumpfile.write("\n")

    return
        
    

def vtkToDict(fname,vtkscalars=['id','type'],
              vtkvectors={'x':('x','y','z'),'F':('Fx','Fy','Fz'),'ux':('ux','uy','uz')},
              integerquantities=['id','type']):
    """
    Convert a vtk file (format vtp usually) to a dictionary which contains arrays
    of data for atoms in the vtk file.
    
    Parameters
    ----------
    fname : string
        Name of the vtk file to read from.
    vtkscalars : list of strings (optional)
        Names of scalar data that should be found in the vtk file. Default is ['id','type'].
    vtkvectors : dictionary (optional)
        Keys are the names of vector data in the vtk file, values are the individual
        component names to be used for the output. Default is 
        {'x':('x','y','z'),'F':('Fx','Fy','Fz'),'ux':('ux','uy','uz')}.
    integerquantities : list of strings (optional)
        Any vtkscalars or values in the vtkvectors dict that are to be treated as integers.
        Default is ['id','type']

    Returns
    -------
    datadict : dictionary 
        Keys will be useful info for e.g. lammps dump output. Values are either numpy arrays,
        ints, or floats.
    natoms : int
        Number of atoms inferred from vtk data.

    """
    
    datadict = {}
    allheaders = []
    vtkheaders = []

    for key,tup in vtkvectors.items():
        for item in tup:
            allheaders.append(item)
        vtkheaders.append(key)

    for item in vtkscalars:
        allheaders.append(item)
        vtkheaders.append(item)

    with open(fname) as vtkfile:
        while allheaders:
            line = vtkfile.readline()

            if line == '':
                raise LookupError (f"Couldn't find the remaining keys {allheaders}")

            for dname in vtkheaders:
                if dname in integerquantities:
                    dtype = 'Int64'
                else:
                    dtype = 'Float64'

                try:
                    ncmp = len(vtkvectors[dname])
                    isvec = True
                except (TypeError, KeyError) as exception:
                    ncmp = 1
                    isvec = False
                    
                matchline = f"<DataArray Name=\"{dname}\" type=\"{dtype}\" NumberOfComponents=\"{ncmp}\" format=\"ascii\">"

                if line.rstrip() == matchline:

                    strarr = vtkfile.readline().rstrip().split()
                    if dtype=='Int64':
                        output = np.array([int(lf) for lf in strarr])
                    else:
                        output = np.array([float(lf) for lf in strarr])
                        
                    if isvec:
                        output = output.reshape(len(output)//ncmp,ncmp)
                        for i,item in enumerate(vtkvectors[dname]):
                            datadict[item] = output[:,i]
                            allheaders.remove(item)
                            

                    else:
                        datadict[dname] = output
                        allheaders.remove(dname)



    for i,(header,val) in enumerate(datadict.items()):
        if i == 0:
            natoms = len(val)
        elif len(val) != natoms:
            raise RuntimeError (f"Number of atoms doesn't match for {key}.")


    return datadict,natoms


if __name__ == "__main__":


    #fname = "swap/dummy_0.vtp"
    #datadict = vtkToDict(fname,vtkscalars=['ID','type','molID'],
    #                         vtkvectors={'x':('x','y','z'),'F':('Fx','Fy','Fz'),'ux':('ux','uy','uz')},
    #                         integerquantities=['ID','type','molID'])


    inputpvdname='swap/dummy.pvd'
    outputdumpname = 'back/dump.lmmpstrj'
    boxbounds = np.array([[-128,128],[-128,128],[-128,128]])
    vtkscalars=['ID','molID','type']
    vtkvectors={'x':('x','y','z'),'F':('Fx','Fy','Fz'),'ux':('ux','uy','uz')}
    integerquantities=['ID','molID','type']

    vtkToDump(inputpvdname,outputdumpname,boxbounds,vtkscalars=vtkscalars,
              vtkvectors=vtkvectors,
              integerquantities=integerquantities,
              headerorder=['ID','molID','type','x','y','z','ux','uy','uz',
                           'Fx','Fy','Fz'])

