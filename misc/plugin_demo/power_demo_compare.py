#!/usr/bin/env python3

import h5py as h5
import numpy as np

fin = h5.File('testsim.out.h5','r')

P_a = fin['/Field/power'][()]
P_b = fin['/Field/plugindemo/my_power'][()]

# To test reporting code
# P_a[1,1]=123

if np.all(P_a==P_b):
    print('objects are identical => OK')
else:
    print('There is at least one different value in the matrices => not OK')

