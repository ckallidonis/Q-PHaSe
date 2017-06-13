#
# lib/Qphase_io.py
#
# Functions regarding I/O
# For reading purposes, a pre-defined format
# is assumed.
#
# Copyright (C) Christos Kallidonis, 2016
# This file is part of Q-PHaSe package
#-------------------------------------------

import sys
import h5py
import numpy as np

sys.path.insert(len(sys.path), '../include/')
from Qphase import *

# Get the hdf5 group path for the two-point function dataset
def get_group2pt(flav,drct,nsrc,msq):
    grp = '/nsrc%02d/msq%04d/%s/%s' % (nsrc,msq,flav,drct)
    return grp
#---------------------------------------------------------


# Get the hdf5 group path for the three-point function dataset
def get_group3pt(op_drct, proj_tag, flav, nsrc, dt, msq):
    grp = '/nsrc%02d/msq%04d/%s/%s/dt%02d/%s' % (nsrc,msq,op_drct,proj_tag,dt,flav)
    return grp
#---------------------------------------------------------

    
# Get the Q-square list for the two-point function
def get_momenta2pt(fname,nsrc,msqMax):
    with h5py.File(fname, 'r') as fp:
        QsqAll = list(fp['nsrc%02d' % nsrc].keys())

    if( msqMax > int(QsqAll[-1].split("msq")[1]) ):
        print('get_momenta2pt: msqMax = %d requested is larger than largest available Qsq = %d. Exiting.' % (msqMax,int(QsqAll[-1].split("msq")[1])) )
        sys.exit()
        
    Qsq = []
    for mom in QsqAll:
        imom = int(mom.split("msq")[1])
        if( imom <= msqMax ):
            Qsq.append(imom)
                
    return Qsq      
#---------------------------------------------------------


# Get the Q-square list for the three-point function
def get_momenta3pt(fname,nsrc,msqMax):
    with h5py.File(fname, 'r') as fp:
        QsqAll = list(fp['nsrc%02d' % nsrc].keys())

    if( msqMax > int(QsqAll[-1].split("msq")[1]) ):
        print('get_momenta3pt: msqMax = %d requested is larger than largest available Qsq = %d. Exiting.' % (msqMax,int(QsqAll[-1].split("msq")[1])) )
        sys.exit()
        
    Qsq = []
    for mom in QsqAll:
        imom = int(mom.split("msq")[1])
        if( imom <= msqMax ):
            Qsq.append(imom)
                
    return Qsq      
#---------------------------------------------------------


# Check if the trajectories for each two-point function case are the same (Qsq, flavor, direction)
def check_traj2pt(twop,traj,Qsq,traj_tmpl):

    traj_tmpl = traj[0][(flav_list2pt[0],drct_list2pt[0])]

    for msq in range(len(Qsq)):
        for flav in flav_list2pt:
            for drct in drct_list2pt:
                if( traj[msq][(flav,drct)] != traj_tmpl ):
                    print('check_traj2pt: Error, trajectories do not agree! Exiting.')
                    sys.exit()
#---------------------------------------------------------


# Check if the trajectories for each three-point function case are the same (Qsq, flavor, operator direction, projector)
def check_traj3pt(thrp,traj,Qsq,basis_tag,optr_tag,proj_list,traj_tmpl):
    
    for msq in range(len(Qsq)):
        for op_drct in operator_tags3pt[optr_tag][basis_tag]:
            for proj_tag in projector_tags3pt[optr_tag][proj_list]:
                for flav in flav_list3pt:
                    if (traj[msq][(op_drct,proj_tag,flav)] != traj_tmpl ):
                        print('check_traj3pt: Error, trajectories do not agree! Exiting.')
                        sys.exit()
#---------------------------------------------------------


# Read the two-point function
def get_twop(fname, nsrc, Qsq):

    twop = []
    mvec = []
    traj = []
        
    with h5py.File(fname, 'r') as fp:
        for msq in range(len(Qsq)):
            twop.append({})
            mvec.append({})
            traj.append({})
            for flav in flav_list2pt:
                for drct in drct_list2pt:
                    grp = get_group2pt(flav,drct,nsrc,Qsq[msq])
                    twop[msq][(flav,drct)] = np.array(fp[grp]['arr'])
                    mvec[msq][(flav,drct)] = np.array(fp[grp]['mvec'])
                    traj[msq][(flav,drct)] = list(fp[grp]['trajs'])

    traj_tmpl = traj[0][(flav_list2pt[0],drct_list2pt[0])] # If all trajectories are the same, just pick the first one
    check_traj2pt(twop,traj,Qsq,traj_tmpl) # Check that the configurations are exactly the same for all flavors, directions and momenta
  
    twopInfo = {}
    twopInfo['Ntraj'] = len(traj_tmpl)
    twopInfo['traj_list'] = traj_tmpl
    twopInfo['T'] = np.shape(twop[0][(flav_list2pt[0],drct_list2pt[0])])[1]
    
    return twop,mvec,traj,twopInfo
#---------------------------------------------------------


# Read the three-point function
def get_thrp(fname, nsrc, tsink, Qsq, basis_tag, optr_tag):

    check_support(basis_tag,optr_tag)
        
    proj_list = get_3ptproj_tag(basis_tag,optr_tag)
    
    thrp = []
    mvec = []
    traj = []    
    with h5py.File(fname, 'r') as fp:
        for msq in range(len(Qsq)):
            thrp.append({})
            mvec.append({})
            traj.append({})
            for op_drct in operator_tags3pt[optr_tag][basis_tag]:
                for proj_tag in projector_tags3pt[optr_tag][proj_list]:
                    for flav in flav_list3pt:
                        grp = get_group3pt(op_drct, proj_tag, flav, nsrc, tsink, Qsq[msq])
                        
                        thrp[msq][(op_drct,proj_tag,flav)] = np.array(fp[grp]['arr'])
                        mvec[msq][(op_drct,proj_tag,flav)] = np.array(fp[grp]['mvec'])
                        traj[msq][(op_drct,proj_tag,flav)] = list(fp[grp]['trajs'])


    tmpl_tpl = (operator_tags3pt[optr_tag][basis_tag][0],projector_tags3pt[optr_tag][proj_list][0],flav_list3pt[0])
    traj_tmpl = traj[0][tmpl_tpl] # If all trajectories are the same, just pick the first one
    check_traj3pt(thrp,traj,Qsq,basis_tag,optr_tag,proj_list,traj_tmpl) # Check that the configurations are exactly the same for all flavors, operator directions, projectors and momenta

  
    thrpInfo = {}
    thrpInfo['Ntraj'] = len(traj_tmpl)
    thrpInfo['traj_list'] = traj_tmpl
    thrpInfo['T'] = np.shape(thrp[0][tmpl_tpl])[1]
    
    return thrp,mvec,traj,thrpInfo
#---------------------------------------------------------
