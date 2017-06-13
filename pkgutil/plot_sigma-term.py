#
# pkgutil/plot_ax-chrg.py
#
# Creates plateau plots for
# the scalar charge
#
# Copyright (C) Christos Kallidonis, 2016
# This file is part of Q-PHaSe package
#----------------------------------------

import sys
import numpy as np
from matplotlib import pyplot as plt

# Read in command line arguments
if ( len(sys.argv) != 6 ):
    print('Usage: %s <file_dir> <flav> <tsinks ..,..,..> <tsink_fit> <fit_range>' % (sys.argv[0]))
    sys.exit()

file_dir = sys.argv[1] 
plt_flv  = sys.argv[2]
ts_lst   = sys.argv[3]
ts_fit   = sys.argv[4]
fit_rng  = sys.argv[5]

obs_tag = 'sigma-term'

ts_tpl = tuple(ts_lst.split(','))


vec_keys = ('val','err','chi')
vep_keys = ('val','err')

Nvec = len(vec_keys)
Nvep = len(vep_keys)

pts_fmt = ('ro','bs','g^')
lne_fmt = ('r-','b-','g-')
cft_col = ('red','blue','green')    


fits = {}
plat = {}
tins = {}
tplt = {}
plat_plt = {}
lgnd_lst = []
for ts in ts_tpl:
    plat_fname = '%s/%s_plat_tsink%s.txt' % (file_dir,obs_tag,ts)
    fits_fname = '%s/%s_fits_tsink%s.txt' % (file_dir,obs_tag,ts)

    # Read the fits data
    with open(fits_fname, 'r') as fp:
        flav_list = fp.readline().strip().split()
        Nflv = len(flav_list)
        
        for line in fp:
            line_lst = line.strip().split()
            rng = line_lst[0]            
            for flav in flav_list:
                fc = flav_list.index(flav)
                for val in vec_keys:
                    vc = vec_keys.index(val)
                    fits[(ts,rng,flav,val)] = float(line_lst[1+vc+Nvec*fc])
                    
    # Read the plateaus
    with open(plat_fname, 'r') as fp:
        flav_list = fp.readline().strip().split()
        Nflv = len(flav_list)

        tins[ts] = []
        tplt[ts] = []
        for flav in flav_list:
            plat[(ts,flav,'val')] = []
            plat[(ts,flav,'err')] = []

        for line in fp:
            line_lst = line.strip().split()
            tins[ts].append(int(line_lst[0]))
            tplt[ts].append(int(line_lst[1]))
            for flav in flav_list:
                fc = flav_list.index(flav)
                for val in vep_keys:
                    vc = vep_keys.index(val)
                    plat[(ts,flav,val)].append(float(line_lst[2+vc+Nvep*fc]))

                    
    tins[ts] = np.array(tins[ts])
    tplt[ts] = np.array(tplt[ts])
    for flav in flav_list:
        for val in vep_keys:
            plat[(ts,flav,val)] = np.array(plat[(ts,flav,val)])

    # Plot the plateaus
    plat_plt[ts] = plt.errorbar(tplt[ts],plat[(ts,plt_flv,'val')],yerr=plat[(ts,plt_flv,'err')], linestyle="None",fmt=pts_fmt[ts_tpl.index(ts)],label=r'$t_s/a = %s$' % (ts))

    # For legend purposes
    lgnd_lst.append(plat_plt[ts])

    
# Plot the constant fit
fit_val = fits[(ts_fit,fit_rng,plt_flv,'val')]
fit_err = fits[(ts_fit,fit_rng,plt_flv,'err')]
fit_vme = fit_val - fit_err
fit_vpe = fit_val + fit_err
tini = int(fit_rng.split('-')[0]) - int(ts_fit)//2
tfin = int(fit_rng.split('-')[1]) - int(ts_fit)//2

cfit  = plt.plot( (tini,tfin), (fit_val,fit_val), lne_fmt[ts_tpl.index(ts_fit)],linewidth=1.5)
cfill = plt.fill_between( (tini,tfin) , fit_vme, fit_vpe, facecolor=cft_col[ts_tpl.index(ts_fit)], alpha=0.35)

# Set the x,y ranges
max_ts = max(ts_tpl)
xmin = min(tplt[max_ts]) - 0.5
xmax = max(tplt[max_ts]) + 0.5
axes = plt.gca()
axes.set_xlim([xmin,xmax])
#axes.set_ylim([10,45])

plt.legend(handles=lgnd_lst,numpoints=1,loc='lower left',ncol=len(ts_tpl))

plt.xlabel(r'$(t_{\rm ins} - t_s/2)/a$')


plt.show()
