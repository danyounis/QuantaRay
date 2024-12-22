# wdb.py
# bin virtual-particle weight data

import h5py
import numpy as np

# load data
f = h5py.File('ed_data.h5','r')

# get & bin virtual-particle weights
wt = np.array(f['data'][5,:])

y, x = np.histogram(wt, bins=int(500), range=[np.min(wt), np.max(wt)])

# output
hf = h5py.File('wdist.h5','w')
hf.create_dataset('y', data=y)
hf.create_dataset('x', data=x)
hf.create_dataset('nde', data=f['data'].shape[-1])
hf.close()
