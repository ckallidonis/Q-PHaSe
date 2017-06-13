#
# axial_charge.py
#
# Reads two- and three-point functions
# and calculate the axial charge.
#
# Copyright (C) Christos Kallidonis, 2016
# This file is part of Q-PHaSe package
#----------------------------------------

import sys
import numpy as np
import h5py
from matplotlib import pyplot as plt

sys.path.insert(len(sys.path), '../lib/')
sys.path.insert(len(sys.path), '../include/')
from gamma_basis import *
from Qphase_io   import *
from Qphase_math import *
from Qphase_jackknife import *
from Qphase_matrix_element import *

# Read in command line arguments
if ( len(sys.argv) != 10 ):
    print('Usage: %s Gamma_basis(QuaHoG, QUDA_UKQCD, QUDA_DeGrand) nsrc msqMax twop-filename thrp-filename tsink Z_renorm binsize output_dir' % (sys.argv[0]))
    sys.exit()
    
basis_tag = sys.argv[1]
nsrc      = sys.argv[2]
msqMax    = sys.argv[3]
fname2pt  = sys.argv[4]
fname3pt  = sys.argv[5]
tsink     = sys.argv[6]
Zfac      = sys.argv[7]
binsize   = sys.argv[8]
outdir    = sys.argv[9]

nsrc    = int(nsrc)
msqMax  = int(msqMax)
tsink   = int(tsink)
Zfac    = float(Zfac)
binsize = int(binsize)

optr_tag = 'axial-vec' # We are computing the axial charge
obs_tag = 'axial_charge'
#--------------------------------------------------------------------------------

g1,g2,g3,g4,g5 = load_basis(basis_tag) # Determine the Gamma-matrices basis to use throughout


Qsq2pt = get_momenta2pt(fname2pt,nsrc,msqMax) # Get the momenta of the two-point function
print('Qsq2pt = ',Qsq2pt)

twop_raw,mvec2pt,traj2pt,twopInfo = get_twop(fname2pt,nsrc,Qsq2pt) # Read the two-point function

twop = average_twop(twop_raw,Qsq2pt,twopInfo,'momenta','ppm','fwd')  # Average two-point function (only over momenta)

print('2pt: Ntraj = %d, T = %d' % (twopInfo['Ntraj'],twopInfo['T']))
#--------------------------------------------------------------------------------

# Get the momenta of the three-point function
Qsq3pt = get_momenta3pt(fname3pt,nsrc,msqMax)
print('Qsq3pt = ',Qsq3pt)

if( Qsq2pt != Qsq3pt ):
    print('axial_charge.py: Q-square lists for two- and three-point functions do not agree! Exiting.')
    sys.exit()
    
# Read the three-point function
thrp_raw,mvec3pt,traj3pt,thrpInfo = get_thrp(fname3pt, nsrc, tsink, Qsq3pt, basis_tag, optr_tag)


# Check that the trajectories of the two- and three-point functions are the same before proceeding
if (twopInfo['traj_list'] != thrpInfo['traj_list']):
    print('axial_charge.py: Trajectory lists for two- and three-point functions do not agree! Exiting.')
    sys.exit()
    
print('3pt: Ntraj = %d, Ntins = %d' % (thrpInfo['Ntraj'],thrpInfo['T']))
#--------------------------------------------------------------------------------

thrp = extractMatrixElement(thrp_raw,thrpInfo,Qsq3pt,basis_tag,optr_tag) # Extract the axial charge from the three-point function


# Perform Jackknife analysis (only on zero-momentum)
msq = Qsq2pt.index(0) # Pick zero-momentum for two-point function
twop_bins, Nbins = jackknife_binning(twop[msq],twopInfo,binsize)
print( '2pt Jackknife: Nbins = %d' % (Nbins) )

thrp_bins = {}    
for comb in flcombs_conn:
    thrp_bins[comb], Nbins = jackknife_binning(thrp[comb],thrpInfo,binsize) 
print( '3pt Jackknife: Nbins = %d' % (Nbins) )


ratio_bin = {}
ratio_ave = {}
T = thrpInfo['T']
for comb in flcombs_conn:
    ratio_bin[comb] = np.zeros((Nbins,T))
    ratio_ave[comb] = {'val' : np.zeros(T) , 'err' : np.zeros(T)}
    for t in range(T):
        for bn in range(Nbins):
            ratio_bin[comb][bn][t] = -thrp_bins[comb][bn][t].imag / twop_bins[bn][tsink].real # Ratio yielding the axial charge

    ratio_bin[comb] *= Zfac # Renormalization!
            
    ratio_ave[comb]['val'] = np.average(ratio_bin[comb],axis=0)                                  # Average ratio over the bins for each time-slice
    ratio_ave[comb]['err'] = jackknife_error(ratio_bin[comb][:,:], ratio_ave[comb]['val'], Nbins, T)
#--------------------------------------------------------------------------------


# Perform constant fits on each bin for all available fit-ranges
T = thrpInfo['T']
nfits  = T//2-1
tini   = 1
tfin   = T-2
rcfit  = {}
chisq  = {}
rcfit_ave = {}
chisq_ave = {}
fit_ranges = {}
for comb in flcombs_conn:
    rcfit[comb] = np.zeros((Nbins,nfits))
    chisq[comb] = np.zeros((Nbins,nfits))
    rcfit_ave[comb] = {'val' : np.zeros(nfits) , 'err' : np.zeros(nfits)}
    chisq_ave[comb] = {'val' : np.zeros(nfits) , 'err' : np.zeros(nfits)}
    fit_ranges[comb] = []
    for nf in range(nfits): # Loop over the various fit ranges
        fit_values = {}
        tstart = tini+nf
        tstop  = tfin-nf
        nPts   = tstart - tstop + 1
        fit_ranges[comb].append('%d-%d' % (tstart,tstop))
        
        for bn in range(Nbins):
            fit_values['val'] = ratio_bin[comb][bn,tstart:tstop+1]
            fit_values['err'] = ratio_ave[comb]['err'][tstart:tstop+1]
            rcfit[comb][bn][nf] = constant_fit(fit_values)                     # Peform the constant fit on each bin
            chisq[comb][bn][nf] = chisq_const(fit_values,rcfit[comb][bn][nf])  # and calculate the chi_square
            
    rcfit_ave[comb]['val'] = np.average(rcfit[comb],axis=0) # Final value of the ratio for each fit-range,      averaged over the fits on the bins
    chisq_ave[comb]['val'] = np.average(chisq[comb],axis=0) # Final value of the chi_square for each fit-range, averaged over the fits on the bins

    rcfit_ave[comb]['err'] = jackknife_error( rcfit[comb][:,:], rcfit_ave[comb]['val'], Nbins, nfits) # Jackknife
    chisq_ave[comb]['err'] = jackknife_error( chisq[comb][:,:], chisq_ave[comb]['val'], Nbins, nfits) # error
#--------------------------------------------------------------------------------


# Write the time-dependence of the ratio
out_fname = '%s/%s_plat_tsink%02d.txt' % (outdir,obs_tag,tsink)

f = open(out_fname,'w')

for comb in flcombs_conn:
    f.write('\t\t%s\t' % (comb))
f.write('\n')

T = thrpInfo['T']
for t in range(T):
    tp = t - tsink/2
    f.write('%2d  %2d\t' % (t,tp))
    for comb in flcombs_conn:
        for v in ('val','err'):
            f.write(' %8.6f' % (ratio_ave[comb][v][t]))
        f.write('\t')
    f.write('\n')

f.close()
#--------------------------------------------------------------------------------


# Write the fit data
out_fname = '%s/%s_fits_tsink%02d.txt' % (outdir,obs_tag,tsink)

f = open(out_fname,'w')

for comb in flcombs_conn:
    f.write('\t\t%s\t\t' % (comb))
f.write('\n')


for fit in range(len(fit_ranges[flcombs_conn[0]])):
    f.write('%s\t' % (fit_ranges[flcombs_conn[0]][fit]))
    for comb in flcombs_conn:
        for v in ('val','err'):
            f.write(' %8.6f' % (rcfit_ave[comb][v][fit]))
        f.write(' %8.6f' % (chisq_ave[comb]['val'][fit]))
        f.write('\t')
    f.write('\n')

f.close()



#################  S C R A T C H  #################

# Testing purposes
# twop_ave = []
# for msq in range(len(Qsq2pt)):
#     twop_ave.append(np.average(twop[msq],axis=0))

# thrp_ave = {}
# for comb in flcombs_conn:
#     thrp_ave[comb] = np.average(thrp[comb],axis=0)


# Ratio = {}
# for comb in flcombs_conn:
#     Ratio[comb] = np.zeros(tsink)
#     for t in range(tsink):
#         Ratio[comb][t] = -thrp_ave[comb][t].imag/twop_ave[0][tsink].real
#         Ratio[comb][t] *= Zfac

