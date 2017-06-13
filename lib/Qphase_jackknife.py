#
# lib/Qphase_jackknife.py
#
# Set of functions for performing
# Jackknife binning
#
# Copyright (C) Christos Kallidonis, 2016
# This file is part of Q-PHaSe package
#----------------------------------------         

import sys
import numpy as np

sys.path.insert(len(sys.path), '../include/')
from Qphase import *


def jackknife_error(bins,ave,Nbins,Npts=1):

    if(Npts==1):
        sqsum = sum(map(lambda x:x*x,ave-bins))
        fac = (Nbins -1) / float(Nbins)
        err = np.sqrt(fac*sqsum)
    else:
        err = np.zeros(Npts)
        for i in range(Npts):
            sqsum = sum(map(lambda x:x*x,ave[i]-bins[:,i]))
            fac = (Nbins -1) / float(Nbins)
            err[i] = np.sqrt(fac*sqsum)

    return err
#---------------------------------------------------------------------------------------


def jackknife_binning(corr, Info, binsize):

    Ntraj = Info['Ntraj']
    T = Info['T']

    # Define the number of bins
    mod   = Ntraj%binsize
    Nbins = int((Ntraj - mod) / binsize)
    bins = np.zeros( (Nbins,T), dtype=np.complex128)

    csum = np.sum(corr,axis=0) # Sum w.r.t to the trajectories
    for t in range(T):
        for m in np.arange(1,mod+1):
            csum[t] -= corr[Ntraj-m][t]  # Throw away data in case there is modulo
            
        for b in range(Nbins):
            bsum = 0
            for k in range(binsize):
                bsum += corr[b*binsize+k][t][0]
                
            bins[b][t] = (csum[t][0] - bsum) / float(Ntraj - binsize - mod) # Bin averages for each bin,t

    return bins, Nbins
#---------------------------------------------------------------------------------------
