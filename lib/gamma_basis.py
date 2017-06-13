#
# lib/gamma_basis.py
#
# Defines the gamma-marices basis
#
# Copyright (C) Christos Kallidonis, 2016
# This file is part of Q-PHaSe package
#----------------------------------------

import sys
from importlib import import_module
sys.path.insert(len(sys.path), '../lib/')
sys.path.insert(len(sys.path), '../include/')
from Qphase import *

def load_basis(basis_tag):

    for key,val in Gamma_bases.items():
        if ( basis_tag == key ):
            print('load_basis: %s basis loaded!' % (basis_tag))
            gammas = import_module("%s" % val)
            break
    else:
        print('load_basis: Usupported basis %s' % (basis_tag))
        sys.exit()
            
    return gammas.g1,gammas.g2,gammas.g3,gammas.g4,gammas.g5
