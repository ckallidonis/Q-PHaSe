#
# lib/Qphase_matrix_element.py
#
# Set of functions for extracting
# matrix elements from three-point functions
#
# Copyright (C) Christos Kallidonis, 2016
# This file is part of Q-PHaSe package
#----------------------------------------

import sys
import numpy as np

sys.path.insert(len(sys.path), '../include/')
from Qphase import *


def extractMatrixElement(thrp_raw,thrpInfo,Qsq,basis_tag,optr_tag):

    check_support(basis_tag,optr_tag)

    curr_lst = operator_tags3pt[optr_tag][basis_tag]
    csgn = optr_sign3pt[basis_tag][optr_tag]

    if(optr_tag == 'axial-vec'):
        msq = Qsq.index(0)
        proj_list = get_3ptproj_tag(basis_tag,optr_tag)
        plist = projector_tags3pt[optr_tag][proj_list]
        thrp = extract_axial_charge(thrp_raw,thrpInfo,curr_lst,plist,csgn) # (momentum-dependence not required here as we are dealing with the axial charge)
    elif(optr_tag == 'scalar'):
        msq = Qsq.index(0)
        proj_list = get_3ptproj_tag(basis_tag,optr_tag)
        plist = projector_tags3pt[optr_tag][proj_list]
        thrp = extract_scalar_charge(thrp_raw,thrpInfo,curr_lst,plist,csgn)
    
        
    return thrp
#---------------------------------------------------------------------------------------


# Function for extracting the scalar charge
def extract_scalar_charge(thrp_raw,thrpInfo,curr_lst,plist,csgn,msq=0):

    thrp = {}
    for flav in flav_list3pt:
        thrp[flav] = thrp_raw[msq][(curr_lst[0],plist[0],flav)]


    thrp['IS'] = thrp['up'] + csgn * thrp['dn']  # Isoscalar combination
    thrp['IV'] = thrp['up'] - csgn * thrp['dn']  # Isovector combination
    
    return thrp
#---------------------------------------------------------------------------------------


# Function for extracting the axial charge
def extract_axial_charge(thrp_raw,thrpInfo,curr_lst,plist,csgn,msq=0):

    pdir  = ['Px','Py','Pz']
    T     = thrpInfo['T']
    Ntraj = thrpInfo['Ntraj']
    NDIR  = 3
    
    # Get the correct projector for each component of the current
    thrp_dir = {}
    thrp = {}
    for flav in flav_list3pt:
        thrp_dir[(curr_lst[0],pdir[0],flav)] = thrp_raw[msq][(curr_lst[0],plist[0],flav)] # x-direction, current g5gx with projector P4 ~ g5gx
        thrp_dir[(curr_lst[1],pdir[1],flav)] = thrp_raw[msq][(curr_lst[1],plist[1],flav)] # y-direction, current g5gy with projector P5 ~ g5gy
        
        # z-direction, for current g5gz, for correct projector, subtract P3-P4-P5 ~ g5(gx+gy+gz) - g5gx - g5gy ~ g5gz
        thrp_dir[(curr_lst[2],pdir[2],flav)] = thrp_raw[msq][(curr_lst[2],plist[2],flav)] - thrp_raw[msq][(curr_lst[2],plist[0],flav)] - thrp_raw[msq][(curr_lst[2],plist[1],flav)]

        thrp[flav] = np.zeros((Ntraj, T, 1), dtype=np.complex128)
        for d in range(NDIR):
            thrp[flav] += thrp_dir[(curr_lst[d],pdir[d],flav)]

        thrp[flav] = thrp[flav] / float(NDIR)

    thrp['IS'] = thrp['up'] + csgn * thrp['dn']  # Isoscalar combination
    thrp['IV'] = thrp['up'] - csgn * thrp['dn']  # Isovector combination
               
    return thrp
#---------------------------------------------------------------------------------------
