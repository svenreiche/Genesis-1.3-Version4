#!/usr/bin/env python3

# Christoph Lechner, Aug 2025
# Demo for automatic testing: reports FEL power at the end of the simulated beamline

import h5py as h5
import sys

# power value at end of simulated beamline (using files in commit id 70be16d; Aug 6, 2025)
Ptest = 0.932785e9

# permissible relative deviation
maxreldev = 0.3

######

f = h5.File('Example1.out.h5','r')
P = f['/Field/power'][()]
f.close()

Pfinal = P[-1]
Pfinal = Pfinal[0]
print(f'FEL power at the end of Example1 run: {Pfinal/1e9:.6f}GW')


# relative deviations exceeding permissible limit result in exit code signalling error
reldev = (Pfinal-Ptest)/Ptest
if abs(reldev)>maxreldev:
    print(f'relative deviation {reldev} exceeds limit of {maxreldev}')
    sys.exit(1)

print(f'OK: relative deviation {reldev} within limit {maxreldev}')
