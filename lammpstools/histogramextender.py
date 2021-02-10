import numpy as np



class extend_order3():
    """
    Given an array type data which has been flattened and
    compressed from shape (N1,N1,N2) to shape
    (N1-2*nskip)*(N1-2*nskip+1)*N2//2, get it back to its
    original shape while putting some filler (to be specified)
    where the values were originally removed.
    
    """

    def __init__(self,Nuv,Na,nskip):

        self.Nuv = Nuv
        self.Na = Na
        self.nskip = nskip

        return

    def flat_length(self):

        tmp = self.Nuv-2*self.nskip
        return tmp*(tmp+1)*self.Na//2

    def flat_index(self,ubin,vbin,abin):
        v_r =  vbin-self.nskip

        u_r = ubin-self.nskip
        N0 = self.Nuv-2*self.nskip+1
        
        return abin+ (v_r + u_r*N0 - (u_r*(u_r+1))//2)*self.Na

        
    def flat_extend(self,ubin,vbin,abin):
        
        return abin + (vbin + ubin*self.Nuv)*self.Na



    def extend(self,flat_array,ncols=1,filler=np.nan):

        if isinstance(filler, (list, tuple, np.ndarray)):
            if len(filler) != ncols:
                raise ValueError("filler must have same length "
                                 "as flat_array[0,:].")
            
        elif isinstance(filler,float):
            # transform filler into array
            filler = [filler]*ncols
            
        
        # shape array to be 2D
        if ncols == 1:
            fl = np.atleast_2d(flat_array).T
        else:
            fl = np.atleast_2d(flat_array)
        
        if fl.shape != (self.flat_length(),ncols):
            raise ValueError("Incorrect input array shape.")
        
        ntotal = self.Nuv*self.Nuv*self.Na

        # array to be returned
        
        extended = np.empty([ntotal,ncols],float)

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


    def reshape(self,flat_array,ncols=1,filler=np.nan):

        extended = self.extend(flat_array,ncols=ncols,filler=filler)

        outs = []
        for col in range(ncols):
            outs.append(extended[:,col].reshape(self.Nuv,self.Nuv,self.Na).T)

        return outs[:ncols]


class extend_order2():
    """
    Given an array type data which has been flattened and
    compressed from shape (N1,N1) to shape
    (N1-2*nskip)*(N1-2*nskip+1)//2, get it back to its
    original shape while putting some filler (to be specified)
    where the values were originally removed.
    
    """

    def __init__(self,Nuv,nskip):

        self.Nuv = Nuv
        self.nskip = nskip

        return

    def flat_length(self):

        tmp = self.Nuv-2*self.nskip
        return tmp*(tmp+1)//2

    def flat_index(self,ubin,vbin):
        v_r =  vbin-self.nskip

        u_r = ubin-self.nskip
        N0 = self.Nuv-2*self.nskip+1
        
        return v_r + u_r*N0 - (u_r*(u_r+1))//2

        
    def flat_extend(self,ubin,vbin):
        
        return vbin + ubin*self.Nuv



    def extend(self,flat_array,ncols=1,filler=np.nan):

        if isinstance(filler, (list, tuple, np.ndarray)):
            if len(filler) != ncols:
                raise ValueError("filler must have same length "
                                 "as flat_array[0,:].")
            
        elif isinstance(filler,float):
            # transform filler into array
            filler = [filler]*ncols
            
        
        # shape array to be 2D
        if ncols == 1:
            fl = np.atleast_2d(flat_array).T
        else:
            fl = np.atleast_2d(flat_array)
        
        if fl.shape != (self.flat_length(),ncols):
            raise ValueError("Incorrect input array shape.")
        
        ntotal = self.Nuv*self.Nuv

        # array to be returned
        
        extended = np.empty([ntotal,ncols],float)

        for col in range(ncols):
            extended[:,col] = filler[col]

        
        for ubin in range(self.nskip,self.Nuv-self.nskip):
            for vbin in range(self.nskip,self.Nuv-ubin):

                    fl = self.flat_index(ubin,vbin)

                    fee = self.flat_extend(ubin,vbin)

                    for col in range(ncols):
                        extended[fee,col] = flat_array[fl,col]



        return extended[:,:ncols]


    def reshape(self,flat_array,ncols=1,filler=np.nan):

        extended = self.extend(flat_array,ncols=ncols,filler=filler)

        outs = []
        for col in range(ncols):
            outs.append(extended[:,col].reshape(self.Nuv,self.Nuv).T)

        return outs[:ncols]

    
    """


    def reduced_flat_index(ubin,vbin,Npos_bins,nskip):
        v_r =  vbin-nskip

        u_r = ubin-nskip
        N0 = Npos_bins-2*nskip+1

        return (v_r + u_r*N0 - (u_r*(u_r+1))//2)

    def reduced_flat_extend(ubin,vbin,Npos_bins,Nangle_bins):

        return (vbin + ubin*Npos_bins)*Nangle_bins

    def v_3D(ulower,uupper,vlower,vupper,alower,aupper):

        v0 = (uupper**3-ulower**3)/3.0
        v0 *= (vupper**3-vlower**3)/3.0
        v0 *= (np.cos(alower) - np.cos(aupper))

        return v0

    def v_2D(ulower,uupper,vlower,vupper,alower,aupper):

        v0 = (uupper**2-ulower**2)/2.0
        v0 *= (vupper**2-vlower**2)/2.0
        v0 *= (aupper-alower)

        return v0


    def compute_threebody(xs,ys,zs,Natoms,npos_bins,nangle_bins,nskip,rc,
                          xL,yL,zL,dimension):
        dr = rc/npos_bins
        dalpha = np.pi/nangle_bins
        hists = np.zeros([nangle_bins,npos_bins,npos_bins],float)                      
        for i in range(Natoms):
            tmpx = xs[i]
            tmpy = ys[i]
            tmpz = zs[i]
            for j in range(Natoms):
                if (i == j):
                    continue
                xij = xs[j]-tmpx
                yij = ys[j]-tmpy
                zij = zs[j]-tmpz

                for k in range(Natoms):
                    if (k == j) or (k==i):
                        continue

                    xik = xs[k]-tmpx
                    yik = ys[k]-tmpy
                    zik = zs[k]-tmpz

                    xjk = xik - xij
                    yjk = yik - yij
                    zjk = zik - zij

                    rij = np.sqrt(xij**2 + yij**2 + zij**2)
                    rik = np.sqrt(xik**2 + yik**2 + zik**2)
                    rjk = np.sqrt(xjk**2 + yjk**2 + zjk**2)


                    ij_bin = int(rij/dr)
                    ik_bin = int(rik/dr)
                    jk_bin = int(rjk/dr)

                    if (ij_bin >= npos_bins or ik_bin >= npos_bins
                        or jk_bin >= npos_bins):
                        continue

                    rij_dot_rik = xij*xik + yij*yik + zij*zik
                    alpha = np.arccos(rij_dot_rik/(rij*rik))



                    alpha_bin = int(alpha/dalpha)
                    hists[alpha_bin,ij_bin,ik_bin] += 1.0


        totbins = (npos_bins-2*nskip)*(npos_bins-2*nskip+1)*nangle_bins//2
        output = np.empty([totbins,Natoms],float)



        if dimension == 2:
            vol = xL*yL
        else:
            vol = xL*yL*zL



        const = vol*vol/(Natoms*Natoms*Natoms)

        for ubin in range(nskip,npos_bins-nskip):
            for vbin in range(nskip,npos_bins-ubin):
                for abin in range(nangle_bins):

                    ulower = ubin*dr
                    uupper = (ubin+1)*dr

                    vlower = vbin*dr
                    vupper = (vbin+1)*dr

                    alower = abin*dalpha
                    aupper = (abin+1)*dalpha


                    if dimension == 2:
                        vf = v_2D(ulower,uupper,vlower,vupper,
                                  alower,aupper)*4*np.pi/const
                    else:
                        vf = v_3D(ulower,uupper,vlower,vupper,
                                  alower,aupper)*8*np.pi**2/const


                    fl = flat_index(ubin,vbin,abin,npos_bins,nangle_bins,nskip)
                    output[fl,0] = int(fl+1)
                    output[fl,1] = (ulower+uupper)/2.0
                    output[fl,2] = (vlower+vupper)/2.0
                    output[fl,3] = (alower+aupper)/2.0
                    output[fl,4] =  hists[abin,vbin,ubin]/vf

        return output


    def extend(output,npos_bins,nangle_bins,nskip):

        ntotal = npos_bins*npos_bins*nangle_bins

        extended = np.zeros([ntotal,5],float)
        for ubin in range(nskip,npos_bins-nskip):
            for vbin in range(nskip,npos_bins-ubin):
                for abin in range(nangle_bins):


                    fl = flat_index(ubin,vbin,abin,npos_bins,nangle_bins,nskip)

                    fee = flat_extend(ubin,vbin,abin,npos_bins,nangle_bins)
                    extended[fee,0] = output[fl,0]
                    extended[fee,1] = output[fl,1]
                    extended[fee,2] = output[fl,2]
                    extended[fee,3] = output[fl,3]
                    extended[fee,4] = output[fl,4]


        return extended

    def int_2d(ext,index,nangle_bins,dalpha):

        i_low = index*nangle_bins
        i_high = (index+1)*nangle_bins

        return np.sum(2*np.cos(ext[i_low:i_high,3])*ext[i_low:i_high,4])*dalpha

    def int_3d(ext,index,nangle_bins,dalpha):

        i_low = index*nangle_bins
        i_high = (index+1)*nangle_bins

        return np.sum(np.cos(ext[i_low:i_high,3])*np.sin(ext[i_low:i_high,3])*ext[i_low:i_high,4])*dalpha

    def contract(extended,npos_bins,nangle_bins,nskip,rc,dimension):

        dr = rc/npos_bins
        dalpha = np.pi/nangle_bins

        totbins = (npos_bins-2*nskip)*(npos_bins-2*nskip+1)//2

        contracted = np.zeros([totbins,4],float)

        if dimension == 2:
            integrand = int_2d
        else:
            integrand = int_3d


        for ubin in range(nskip,npos_bins-nskip):
            for vbin in range(nskip,npos_bins-ubin):
                index = vbin+npos_bins*ubin
                tmpout = integrand(extended,index,nangle_bins,dalpha)
                fl = reduced_flat_index(ubin,vbin,npos_bins,nskip)

                fee = reduced_flat_extend(ubin,vbin,npos_bins,nangle_bins)
                contracted[fl,0] = fl+1
                contracted[fl,1] = extended[fee,1]
                contracted[fl,2] = extended[fee,2]
                contracted[fl,3] = tmpout


        return contracted
"""


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


    
    eo3 = extend_order3(Nuv,Na,nskip)

    dum = np.empty([eo3.flat_length(),3],float)


    for ubin in range(nskip,Nuv-nskip):
        for vbin in range(nskip,Nuv-ubin):
            for abin in range(Na):

                fl = eo3.flat_index(ubin,vbin,abin)

                dum[fl,0] = us[ubin]
                dum[fl,1] = vs[vbin]
                dum[fl,2] = angles[abin]


    full = eo3.reshape(dum,ncols=3)


    print([full[0][0,:,:],full[1][0,:,:]])

    print(np.meshgrid(us,vs))



    # check order 2 tensor
    eo2 = extend_order2(Nuv,nskip)

    dum = np.empty([eo2.flat_length(),2],float)


    for ubin in range(nskip,Nuv-nskip):
        for vbin in range(nskip,Nuv-ubin):

            fl = eo2.flat_index(ubin,vbin)
            
            dum[fl,0] = us[ubin]
            dum[fl,1] = vs[vbin]


    full = eo2.reshape(dum,ncols=2)


    print(full)


    print(np.meshgrid(us,vs))

