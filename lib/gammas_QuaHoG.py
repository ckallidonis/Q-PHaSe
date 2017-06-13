#
# lib/gammas_QuaHoG.py
#
# Contains definitions of gamma-matrices
# according to the QuaHoG (same as tmLQCD) chiral basis.
#
# Copyright (C) Christos Kallidonis, 2016
# This file is part of Q-PHaSe package
#----------------------------------------

import numpy as np

I = 0 + 1j  # Imaginary unit

g0 = np.matrix([[ 0, 0,-1, 0 ],
                [ 0, 0, 0,-1 ],
                [-1, 0, 0, 0 ],
                [ 0,-1, 0, 0 ]])

g1 = np.matrix([[ 0, 0, 0,-I ],
                [ 0, 0,-I, 0 ],
                [ 0, I, 0, 0 ],
                [ I, 0, 0, 0 ]])

g2 = np.matrix([[ 0, 0, 0,-1 ],
                [ 0, 0, 1, 0 ],
                [ 0, 1, 0, 0 ],
                [-1, 0, 0, 0 ]])

g3 = np.matrix([[ 0, 0,-I, 0 ],
                [ 0, 0, 0, I ],
                [ I, 0, 0, 0 ],
                [ 0,-I, 0, 0 ]])

g4 = np.matrix([[ 0, 0, 1, 0 ],
                [ 0, 0, 0, 1 ],
                [ 1, 0, 0, 0 ],
                [ 0, 1, 0, 0 ]])

g5 = np.matrix(g1 * g2 * g3 * g4)
