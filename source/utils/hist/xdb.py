# xdb.py
# bin virtual-particle position data

import h5py
import numpy as np

# distribution parameters
bins = int(1e+3)
prange = [[-2e+3, 2e+3], [-2e+3, 2e+3]] # au

# load data
f = h5py.File('ed_data.h5','r')

# create container variables
x  = np.array(f['data'][0,:]) # data dims: (6,nde)
y  = np.array(f['data'][1,:])
ph = np.array(f['data'][4,:])
wt = np.array(f['data'][5,:])

# bin particles
weights = wt
W0, xedges, yedges = np.histogram2d(x, y, bins=bins, range=prange, weights=weights)

weights = np.sqrt(wt)*np.cos(ph)
ReW, xedges, yedges = np.histogram2d(x, y, bins=bins, range=prange, weights=weights)

weights = np.sqrt(wt)*np.sin(ph)
ImW, xedges, yedges = np.histogram2d(x, y, bins=bins, range=prange, weights=weights)

# output
hf = h5py.File('xdist.h5','w')
hf.create_dataset('W0', data=W0)
hf.create_dataset('Re(W)', data=ReW)
hf.create_dataset('Im(W)', data=ImW)
hf.create_dataset('x1', data=xedges)
hf.create_dataset('x2', data=yedges)
hf.create_dataset('nde', data=f['data'].shape[-1])
hf.close()
