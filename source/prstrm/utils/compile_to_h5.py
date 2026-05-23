# compile_to_h5.py
import h5py
import numpy as np
from pathlib import Path

mpath = Path('../')
output_fname = 'output.h5'

# save-data precision
precx = {}
precx['float32'] = np.int8
precx['float64'] = np.float32
precx['int32'] = np.int8

out = {}
out['lepton/charge'] = np.fromfile(Path.joinpath(mpath,'charge_lepton.dat'), dtype=np.float32).astype(precx['float32'])
out['lepton/energy'] = np.fromfile(Path.joinpath(mpath,'energy_lepton.dat'), dtype=np.float64).astype(precx['float64'])
out['lepton/quantum parameter'] = np.fromfile(Path.joinpath(mpath,'qparam_lepton.dat'), dtype=np.float64).astype(precx['float64'])
out['lepton/creation time'] = np.fromfile(Path.joinpath(mpath,'x0_lepton.dat'), dtype=np.float64).astype(precx['float64'])
out['lepton/creation energy'] = np.fromfile(Path.joinpath(mpath,'p0_lepton.dat'), dtype=np.float64).astype(precx['float64'])

x1f = np.fromfile(Path.joinpath(mpath,'x1f_lepton.dat'), dtype=np.float64).astype(precx['float64'])
x2f = np.fromfile(Path.joinpath(mpath,'x2f_lepton.dat'), dtype=np.float64).astype(precx['float64'])
x3f = np.fromfile(Path.joinpath(mpath,'x3f_lepton.dat'), dtype=np.float64).astype(precx['float64'])

p1f = np.fromfile(Path.joinpath(mpath,'p1f_lepton.dat'), dtype=np.float64).astype(precx['float64'])
p2f = np.fromfile(Path.joinpath(mpath,'p2f_lepton.dat'), dtype=np.float64).astype(precx['float64'])
p3f = np.fromfile(Path.joinpath(mpath,'p3f_lepton.dat'), dtype=np.float64).astype(precx['float64'])

out['lepton/final position'] = np.stack((x1f,x2f,x3f), axis=0)
out['lepton/final momentum'] = np.stack((p1f,p2f,p3f), axis=0)

out['photon/energy'] = np.fromfile(Path.joinpath(mpath,'energy_photon.dat'), dtype=np.float64).astype(precx['float64'])
out['photon/quantum parameter'] = np.fromfile(Path.joinpath(mpath,'qparam_photon.dat'), dtype=np.float64).astype(precx['float64'])
out['photon/creation time'] = np.fromfile(Path.joinpath(mpath,'x0_photon.dat'), dtype=np.float64).astype(precx['float64'])

out['annihilated photon/energy'] = np.fromfile(Path.joinpath(mpath,'energy_annihilated-photon.dat'), dtype=np.float64).astype(precx['float64'])
out['annihilated photon/quantum parameter'] = np.fromfile(Path.joinpath(mpath,'qparam_annihilated-photon.dat'), dtype=np.float64).astype(precx['float64'])
out['annihilated photon/creation time'] = np.fromfile(Path.joinpath(mpath,'x0_annihilated-photon.dat'), dtype=np.float64).astype(precx['float64'])

x1f = np.fromfile(Path.joinpath(mpath,'x1f_annihilated-photon.dat'), dtype=np.float64).astype(precx['float64'])
x2f = np.fromfile(Path.joinpath(mpath,'x2f_annihilated-photon.dat'), dtype=np.float64).astype(precx['float64'])
x3f = np.fromfile(Path.joinpath(mpath,'x3f_annihilated-photon.dat'), dtype=np.float64).astype(precx['float64'])

p1f = np.fromfile(Path.joinpath(mpath,'p1f_annihilated-photon.dat'), dtype=np.float64).astype(precx['float64'])
p2f = np.fromfile(Path.joinpath(mpath,'p2f_annihilated-photon.dat'), dtype=np.float64).astype(precx['float64'])
p3f = np.fromfile(Path.joinpath(mpath,'p3f_annihilated-photon.dat'), dtype=np.float64).astype(precx['float64'])

out['annihilated photon/final position'] = np.stack((x1f,x2f,x3f), axis=0)
out['annihilated photon/final momentum'] = np.stack((p1f,p2f,p3f), axis=0)

out['lepton/id'] = np.fromfile(Path.joinpath(mpath,'id_lepton.dat'), dtype=np.int32).astype(precx['int32'])
out['photon/id'] = np.fromfile(Path.joinpath(mpath,'id_photon.dat'), dtype=np.int32).astype(precx['int32'])
out['annihilated photon/id'] = np.fromfile(Path.joinpath(mpath,'id_annihilated-photon.dat'), dtype=np.int32).astype(precx['int32'])

with h5py.File(Path.joinpath(mpath,output_fname),'w') as hf:
    for key in out:
        hf.create_dataset(key, data=out[key])
