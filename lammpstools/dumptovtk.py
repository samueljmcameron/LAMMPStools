from .dumploader import DumpLoader
import numpy as np



def DumpToVTK(filein,fileout,stepmax = 10,stepmin=0,
              dt=1e-5,integerquantities=['id','type'],
              vtkscalars=['id','type'],
              vtkvectors={'F':('Fx','Fy','Fz'),'ux':('ux','uy','uz')},
              lowestid=0,idlocation=0):
    """
    Function to convert dump file from lammps output to a series of vtkfiles
    and a collection to track the timesteps. 

    Parameters
    ----------
    filein : string
        Name of the dump file to be loaded in.
    fileout : string
        Desired names of the vtkoutput files, without any suffix (e.g. excluded .vtp)
    stepmax : int
        Maximum step to output to vtk files. Default is 10.
    stepmax : int
        Minimum step to output to vtk files. Default is 0.
    integerquantities : list of strings
        Names of all quantities that should be stored as integers
        instead of floats. Default is ['id','type'].
    vtkscalars : list of strings
        Names of all scalar quantities to be saved to vtk files.
        Default is ['id','type'].
    vtkvectors: dictionary with strings for keys and tuple of strings for values
        Names of all vector quantities to be saved to vtk files (keys)
        and a tuple of the names of the components of these 
        quantities to be found in the dump file. Default is
        {'F' : ('Fx','Fy','Fz'), 'ux' : ('ux','uy','uz')}
    lowestid : int (optional)
        Lowest id that an atom can have. Default is 1.
    idlocation: int >= 0 (optional)
        In the header list of the dump file, the location of the keyword that labels
        the atoms ('id'). Default is zero, meaning the first word in the header is 'id'
        or similar.

    Returns
    -------
    Nothing, just saves files.

    """

    if (stepmin > stepmax):
        raise ValueError("Cannot have minimum step larger than maximum step.")
    dump = DumpLoader(filein,integerquantities=integerquantities,lowestid=lowestid,
                      idlocation=idlocation,min_step=stepmin,max_step=stepmax)

    timesteps=[]
    fnames=[]
    
    for datastep in dump.data:

        if (datastep['step'] > stepmax):
            break
        if (datastep['step'] < stepmin):
            continue

        numpoints = datastep['N']
        timesteps.append(datastep['step'])
        fnames.append(fileout + "_" + str(timesteps[-1]) + ".vtp")
        fname = fnames[-1]

        header = "<?xml version=\"1.0\"?>\n"
        header += "<VTKFile type=\"PolyData\"  version=\"1.0\" byte_order=\"LittleEndian\">\n"
        header += "<PolyData>\n"
        header += "<Piece NumberOfPoints=\"" + str(numpoints)
        header += "\" NumberOfLines=\"1\">\n"
        
        with open(fname,"w") as fout:
            fout.write(header)
            fout.write("<Points>\n")
            xs = datastep['x']
            ys = datastep['y']
            zs = datastep['z']
            positions = np.vstack((xs,ys,zs)).transpose()
            fout.write("<DataArray Name=\"x\" type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n") 
            np.savetxt(fout,positions.flatten(),newline=' ',fmt='%f')
            fout.write("\n</DataArray>\n")
            fout.write("</Points>\n")
            fout.write("<PointData>\n")


            for vname,tup in vtkvectors.items():
                
                vec = ()
                item_type = "Float64"
                for key in tup:
                    vec = (*vec,datastep[key])
                    if key in integerquantities:
                        item_type = "Int64"
                
                vec = np.vstack(vec).transpose()
                fout.write(f"<DataArray Name=\"{vname}\" type=\"{item_type}\" NumberOfComponents=\"{len(tup)}\" format=\"ascii\">\n")
                if item_type == "Int64":
                    fmt = '%d'
                else:
                    fmt = '%f'
                np.savetxt(fout,vec.flatten(),newline=' ',fmt=fmt)
                fout.write("\n</DataArray>\n")


            for key in vtkscalars:
                item_type = "Float64"
                if key in integerquantities:
                    item_type = "Int64"
                scal = datastep[key]
                if item_type == "Int64":
                    fmt = '%d'
                else:
                    fmt = '%f'                
                fout.write(f"<DataArray Name=\"{key}\" type=\"{item_type}\" NumberOfComponents=\"1\" format=\"ascii\">\n")
                np.savetxt(fout,scal,newline=' ',fmt=fmt)
                fout.write("\n</DataArray>\n")


            fout.write("</PointData>\n</Piece>\n</PolyData>\n</VTKFile>")



    header = "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\"  version=\"1.0\""
    header += " byte_order=\"LittleEndian\">\n<Collection>\n"

    with open(fileout+".pvd","w") as fout:

        fout.write(header)

        for tstep,fname in zip(timesteps,fnames):

            nodirfname = fname[1+fname.rfind("/"):]
            if nodirfname == "":
                nodirfname = fname

            line = f"<DataSet timestep=\"{tstep*dt}\" group=\"\" part=\"0\""
            line += f" file=\"{nodirfname}\"/>\n"
            fout.write(line)

        fout.write("</Collection>\n</VTKFile>")



    return;
