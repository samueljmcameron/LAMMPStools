import numpy as np



class CondensedArray_oThree():
    """
    Methods to manipulate an array type data which has been
    flattened and compressed from shape (Na,Nuv,Nuv) to a 1D
    array of length (Nuv-2*nskip)*(Nuv-2*nskip+1)*Na//2.

    Attributes
    ----------
    Nuv : int
        Length of final two axes of original array.
    Na : int
        Length of first axis of original array.
    nskip : int
        Number of leading entries which were skipped
        in the final two axes or original array.


    Methods
    -------

    flat_length(self)

        Compute the expected length of the condensed array
        given the attributes of the class.

    flat_index(self,ubin,vbin,abin)

        Compute the flattened index of the condensed array,
        given the indices of the original array, i.e. a one-way
        map from original[abin,vbin,ubin] to
        condensed[flat_index(ubin,vbin,abin)].


        
    flat_extend(self,ubin,vbin,abin)

        Compute the flattened index of the original
        (uncondensed) array, given the indices of this array.


    extend(self,flat_array,ncols=1,filler=np.nan)

        Extend the condensed, flattened array to the
        uncondensed, flattened array, filling in the
        missing values with a specified value (the
        filler parameter).


    reshape(self,flat_array,ncols=1,filler=np.nan)    

        Reshape the condensed, flattened array to 
        uncondensed, full order two arrays, filling in the
        missing values with a specified value (the
        filler parameter).

    """

    def __init__(self,Nuv,Na,nskip):

        """
        Initialise attributes.

        Parameters
        ----------
        Nuv : int
            Length of final two axes of original array.
        Na : int
            Length of first axis of original array.
        nskip : int
            Number of leading entries which were skipped
            in the final two axes of original array.

        """

        self.Nuv = Nuv
        self.Na = Na
        self.nskip = nskip

        return

    def flat_length(self):

        """
        Compute the expected length of the condensed array
        given the attributes of the class.

        Parameters
        ----------
        None

        Returns
        -------
        out : int
            expected length of the condensed array, i.e.

                (Nuv-2*nskip)*(Nuv-2*nskip+1)*Na//2        

        """
        
        tmp = self.Nuv-2*self.nskip
        return tmp*(tmp+1)*self.Na//2

    def flat_index(self,ubin,vbin,abin):

        """

        Compute the flattened index of the condensed array,
        given the indices of the original array, i.e. a one-way
        map from original[abin,vbin,ubin] to
        condensed[flat_index(ubin,vbin,abin)].

        Parameters
        ----------
        ubin : int
            Index of third axis in original array (takes same
            range of values as index of second axes, vbin).
        vbin : int
            Index of second axis in original array (takes same
            range of values as index of third axes, ubin).
        abin :int
            Index of first axis in original array (does not
            necessarily take same range of values as ubin and
            vbin).
        
        Returns
        -------
        index : int
            Index of flattened, condensed array (typically array 
            of three body data).

        """

        v_r =  vbin-self.nskip

        u_r = ubin-self.nskip
        N0 = self.Nuv-2*self.nskip+1
        
        return abin+ (v_r + u_r*N0 - (u_r*(u_r+1))//2)*self.Na

        
    def flat_extend(self,ubin,vbin,abin):

        """

        Compute the flattened index of the original
        (uncondensed) array, given the indices of this array.

        Parameters
        ----------
        ubin : int
            Index of third axis in original array (takes same
            range of values as index of second axes, vbin).
        vbin : int
            Index of second axis in original array (takes same
            range of values as index of third axes, ubin).
        abin :int
            Index of first axis in original array (does not
            necessarily take same range of values as ubin and
            vbin).
        
        Returns
        -------
        index : int
            Index of flattened, original array .

        """
        
        
        return abin + (vbin + ubin*self.Nuv)*self.Na



    def _get_ncols(self,flat_array):
        
        """
        Get the number of columns in flat_array. 

        Parameters
        ----------
        flat_array : np.array of dimension 1 or 2.
            Error will be thrown in flat_array.ndim
            does not return 1 or 2.

        Returns
        -------
        ncols : int
            Number of columns in flat_array (returns 1
            if flat_array.ndim == 1).
        
        """
        
        if flat_array.ndim == 1:
            ncols = 1

        elif flat_array.ndim == 2:
            ncols = flat_array.shape[1]
        else:
            raise ValueError('flat_array must be 1D or 2D.')

        return ncols
    
    def extend(self,flat_array,filler=np.nan):

        """

        Extend the condensed, flattened array to the
        uncondensed, flattened array, filling in the
        missing values with a specified value (the
        filler parameter).

        Parameters
        ----------
        flat_array : np.ndarray
            The condensed data that was e.g. output by simulation.
            Must be either 1D or 2D, with condensed data along
            the first axis.
        filler : float or array-like 
            Value to fill in for the elements of the full array
            that were omitted in the condensed array. ncols is
            the number of columns in flat_array. Default
            is np.nan.
        
        Returns
        -------
        extended : np.ndarray
            Flattened, full array of shape
            (self.Nuv*self.Nuv*self.Na,) (for 1D input
            flat_array) or, (self.Nuv*self.Nuv*self.Na,k)
            where k = flat_array.shape[1] (for 2D input array).

        """

        ncols = self._get_ncols(flat_array)

        if ncols == 1:
            flat_array = np.atleast_2d(flat_array).T
            
        if isinstance(filler, (list, tuple, np.ndarray)):
            if len(filler) != ncols:
                raise ValueError("filler must have same length "
                                 "as flat_array[0,:].")
            
        elif isinstance(filler,float):
            # transform filler into array
            filler = [filler]*ncols
            
        
        
        if flat_array.shape != (self.flat_length(),ncols):
            raise ValueError("Incorrect input array shape.")
        
        ntotal = self.Nuv*self.Nuv*self.Na

        # array to be returned
        
        extended = np.empty([ntotal,ncols],float)

        # fill in array with filler values.
        for col in range(ncols):
            extended[:,col] = filler[col]

        
        for ubin in range(self.nskip,self.Nuv-self.nskip):
            for vbin in range(self.nskip,self.Nuv-ubin):
                for abin in range(self.Na):

                    fl = self.flat_index(ubin,vbin,abin)

                    fee = self.flat_extend(ubin,vbin,abin)

                    for col in range(ncols):
                        extended[fee,col] = flat_array[fl,col]



        return extended[:,:ncols]


    def reshape(self,flat_array,filler=np.nan):

        """

        Reshape the condensed, flattened array to 
        uncondensed, full order two arrays, filling in the
        missing values with a specified value (the
        filler parameter).

        Parameters
        ----------
        flat_array : np.ndarray
            The condensed data that was e.g. output by simulation.
            Must be either 1D or 2D, with condensed data along
            the first axis.
        filler : float or array-like of length ncols
            Value to fill in for the elements of the full array
            that were omitted in the condensed array. ncols is
            the number of columns in flat_array. Default
            is np.nan.
        
        Returns
        -------
        out : 3D np.array or list (of length ncols) of 3D np.arrays 
            If flat_array is 1D, then return a single 3D array of
            shape (self.Na,self.Nuv,self.Nuv). If flat_array is 2D,
            return list with each item in the list being a 3D array
            of this shape. ncols = flat_array.shape[1] (or 1 if
            the former doesn't exist).

        Examples
        --------
        If the parameter flat_array was a condensed version
        of several 200 x 400 x 400 arrays, then this function would
        return the original arrays as items in a list, but with all
        values that were omitted from flat_array taking the value of
        the parameter filler.

        """        

        ncols = self._get_ncols(flat_array)
        
        extended = self.extend(flat_array,filler=filler)
        
        if ncols == 1:
            outs = extended.reshape(self.Nuv,self.Nuv).T
        else:
            outs = []
            for col in range(ncols):
                outs.append(extended[:,col].reshape(self.Nuv,self.Nuv,self.Na).T)
                
        return outs


class CondensedArray_oTwo(CondensedArray_oThree):

    """

    Methods to manipulate an array type data which has been
    flattened and compressed from shape (Nuv,Nuv) to a 1D
    array of length (Nuv-2*nskip)*(Nuv-2*nskip+1)//2.

    Attributes
    ----------
    Nuv : int
        Length of two axes of original array.
    nskip : int
        Number of leading entries which were skipped
        in the two axes of original array.


    Methods
    -------

    flat_index(self,ubin,vbin)

        Compute the flattened index of the condensed array,
        given the indices of the original array, i.e. a one-way
        map from original[vbin,ubin] to
        condensed[flat_index(ubin,vbin)].

    flat_extend(self,ubin,vbin)

        Compute the flattened index of the original
        (uncondensed) array, given the indices of this array.

    extend(self,flat_array,filler=np.nan)

        Extend the condensed, flattened array to the
        uncondensed, flattened array, filling in the
        missing values with a specified value (the
        filler parameter).

    reshape(self,flat_array,filler=np.nan)

        Reshape the condensed, flattened array to 
        uncondensed, full order two arrays, filling in the
        missing values with a specified value (the
        filler parameter).

    as well as functions inherited from parent class
    CondensedArray_oThree.
    
    """

    def __init__(self,Nuv,nskip):

        """
        Initialise attributes.

        Parameters
        ----------
        Nuv : int
            Length of two axes of original array.
        nskip : int
            Number of leading entries which were skipped
            in the two axes of original array.

        """

        super().__init__(Nuv,1,nskip)
        return


    def flat_index(self,ubin,vbin):

        """

        Compute the flattened index of the condensed array,
        given the indices of the original array, i.e. a one-way
        map from original[vbin,ubin] to
        condensed[flat_index(ubin,vbin)].

        Parameters
        ----------
        ubin : int
            Index of third axis in original array (takes same
            range of values as index of second axes, vbin).
        vbin : int
            Index of second axis in original array (takes same
            range of values as index of third axes, ubin).
        
        Returns
        -------
        index : int
            Index of flattened, condensed array (typically array 
            of three body data).

        """

        return super().flat_index(ubin,vbin,0)
        
    def flat_extend(self,ubin,vbin):

        """

        Compute the flattened index of the original
        (uncondensed) array, given the indices of this array.

        Parameters
        ----------
        ubin : int
            Index of third axis in original array (takes same
            range of values as index of second axes, vbin).
        vbin : int
            Index of second axis in original array (takes same
            range of values as index of third axes, ubin).
        
        Returns
        -------
        index : int
            Index of flattened, original array .

        """


        return super().flat_extend(ubin,vbin,0)

    def extend(self,flat_array,filler=np.nan):

        """

        Extend the condensed, flattened array to the
        uncondensed, flattened array, filling in the
        missing values with a specified value (the
        filler parameter).

        Parameters
        ----------
        flat_array : np.ndarray
            The condensed data that was e.g. output by simulation.
            Must be either 1D or 2D, with condensed data along
            the first axis.
        filler : float or array-like 
            Value to fill in for the elements of the full array
            that were omitted in the condensed array. ncols is
            the number of columns in flat_array. Default
            is np.nan.
        
        Returns
        -------
        extended : np.ndarray
            Flattened, full array of shape
            (self.Nuv*self.Nuv,) (for 1D input
            flat_array) or, (self.Nuv*self.Nuv,k)
            where k = flat_array.shape[1] (for 2D input array).

        """

        ncols = self._get_ncols(flat_array)

        if ncols == 1:
            flat_array = np.atleast_2d(flat_array).T
            
        if isinstance(filler, (list, tuple, np.ndarray)):
            if len(filler) != ncols:
                raise ValueError("filler must have same length "
                                 "as flat_array[0,:].")
            
        elif isinstance(filler,float):
            # transform filler into array
            filler = [filler]*ncols
            
        
        
        if flat_array.shape != (self.flat_length(),ncols):
            raise ValueError("Incorrect input array shape.")
        
        ntotal = self.Nuv*self.Nuv*self.Na

        # array to be returned
        
        extended = np.empty([ntotal,ncols],float)

        # fill in array with filler values.
        for col in range(ncols):
            extended[:,col] = filler[col]

        
        for ubin in range(self.nskip,self.Nuv-self.nskip):
            for vbin in range(self.nskip,self.Nuv-ubin):

                fl = self.flat_index(ubin,vbin)

                fee = self.flat_extend(ubin,vbin)

                for col in range(ncols):
                    extended[fee,col] = flat_array[fl,col]



        return extended[:,:ncols]


    def reshape(self,flat_array,filler=np.nan):

        """

        Reshape the condensed, flattened array to 
        uncondensed, full order two arrays, filling in the
        missing values with a specified value (the
        filler parameter).

        Parameters
        ----------
        flat_array : np.ndarray
            The condensed data that was e.g. output by simulation.
            Must be either 1D or 2D, with condensed data along
            the first axis.
        filler : float or array-like of length ncols
            Value to fill in for the elements of the full array
            that were omitted in the condensed array. ncols is
            the number of columns in flat_array. Default
            is np.nan.
        
        Returns
        -------
        out : 3D np.array or list (of length ncols) of 3D np.arrays 
            If flat_array is 1D, then return a single 3D array of
            shape (self.Nuv,self.Nuv). If flat_array is 2D,
            return list with each item in the list being a 3D array
            of this shape. ncols = flat_array.shape[1] (or 1 if
            the former doesn't exist).

        Examples
        --------
        If the parameter flat_array was a condensed version
        of several 400 x 400 arrays, then this function would
        return the original arrays as items in a list, but with all
        values that were omitted from flat_array taking the value of
        the parameter filler.

        """        

        ncols = self._get_ncols(flat_array)
        
        extended = self.extend(flat_array,filler=filler)

        if ncols == 1:
            outs = extended.reshape(self.Nuv,self.Nuv).T
        else:
            outs = []
            for col in range(ncols):
                outs.append(extended[:,col].reshape(self.Nuv,self.Nuv).T)

        return outs

    


if __name__ == "__main__":

    # check order 3 tensor extender
    Nuv = 4
    Na = 2
    nskip = 0

    rc = 16.0

    dr = rc/Nuv

    da = np.pi/Na


    us = np.linspace(0,rc,num=Nuv,endpoint=False)+0.5*dr
    vs = np.linspace(0,rc,num=Nuv,endpoint=False)+0.5*dr

    angles = np.linspace(0,np.pi,num=Na,endpoint=False)+0.5*da


    
    eo3 = CondensedArray_oThree(Nuv,Na,nskip)

    dum = np.empty([eo3.flat_length(),3],float)


    for ubin in range(nskip,Nuv-nskip):
        for vbin in range(nskip,Nuv-ubin):
            for abin in range(Na):

                fl = eo3.flat_index(ubin,vbin,abin)

                dum[fl,0] = us[ubin]
                dum[fl,1] = vs[vbin]
                dum[fl,2] = angles[abin]


    full = eo3.reshape(dum)


    print([full[0][0,:,:],full[1][0,:,:]])

    print(np.meshgrid(us,vs))



    # check order 2 tensor
    eo2 = CondensedArray_oTwo(Nuv,nskip)

    dum = np.empty([eo2.flat_length(),2],float)


    for ubin in range(nskip,Nuv-nskip):
        for vbin in range(nskip,Nuv-ubin):

            fl = eo2.flat_index(ubin,vbin)
            
            dum[fl,0] = us[ubin]
            dum[fl,1] = vs[vbin]


    full = eo2.reshape(dum)


    print(full)


    print(np.meshgrid(us,vs))

