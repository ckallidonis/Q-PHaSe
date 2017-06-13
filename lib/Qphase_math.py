#
# lib/Qphase_math.py
#
# Set of functions for performing
# mathematical operations
#
# Copyright (C) Christos Kallidonis, 2016
# This file is part of Q-PHaSe package
#----------------------------------------

import sys
import numpy as np

sys.path.insert(len(sys.path), '../include/')
from Qphase import *


# Average the two-point function over the momenta for each Qsq,
# flavor (ppm,pmm) and direction (fwd,bwd)
def average_twop(twop_raw,Qsq,twopInfo,avg_type="all",pflv='',pdrc=''):

    Ntraj = twopInfo['Ntraj']
    T     = twopInfo['T']
    nQsq  = len(Qsq)
    
    twop = np.zeros((nQsq, Ntraj, T, 1),dtype=np.complex128)   
    
    if(avg_type=="all"):
        twop_drct = []
        twop_drct_aver = []
        for msq in range(nQsq):
            twop_drct.append({})
            twop_drct_aver.append({})
            for drct in drct_list2pt:
                twop_drct[msq][drct] = []
                for flav in flav_list2pt:
                    twop_drct[msq][drct].append( np.average(twop_raw[msq][(flav,drct)],axis=2) ) # Average over the momenta for each Qsq, for each flavor and direction
        
                twop_drct_aver[msq][drct] = np.average(twop_drct[msq][drct],axis=0)  # Average over the flavors, for each Qsq and direction

            # Average over forward, backward directions
            for ic in range(Ntraj):
                for it in range(T):
                    twop[msq][ic][it][0] = 0.5 * ( twop_drct_aver[msq]['fwd'][ic][it] - twop_drct_aver[msq]['bwd'][ic][(T-it)%T] )

    elif(avg_type=="momenta"):
        if(pflv=='' or pdrc==''):
            print( 'average_twop: Variables pflv and pdrc must be defined when avg_type = "%s"' % (avg_type) )
            sys.exit()

        print( 'average_twop: Averaging over momenta, flavor = %s , direction = %s' % (pflv,pdrc) )
        for msq in range(nQsq):
            twop[msq] = np.reshape( np.average(twop_raw[msq][(pflv,pdrc)],axis=2)  , (Ntraj,T,1) ) # Average over the momenta for each Qsq, for flavor 'ppm' and direction 'fwd'

            
    return twop
#---------------------------------------------------------------------------------------
    
    
# Function which performs a constant fit
def constant_fit(values):

    Sy = sum( map(lambda x:x,values['val']/(values['err']**2)) )  # Sy = SUM_i val[i]/err[i]**2
    S  = sum( map(lambda x:1.0/(x*x),values['err']) )             # S  = SUM_i ( 1/err[i] )^2

    return Sy/S
#---------------------------------------------------------------------------------------


# Function which calculates the chi_square for a constant fit
def chisq_const(values, fit):

    return sum( map(lambda x:x,((values['val']-fit)/values['err'])**2) )
#---------------------------------------------------------------------------------------
