# dat2h5.py
# Convert Fortran raw binary data to HDF5 format.
#
# Since CUDA Fortran doesn't support HDF5, this script takes the .dat files from
# QuantaRay's CUDA QCR and packages them into the standard QM H5 file.
#
# An example of an unformatted sequential file in Fortran would be written as:
#   OPEN(1, FILE=filename, FORM='unformatted')
#   WRITE(1) variable
#
# Note: Quantities are transposed to match the dimensionality swapping between
# Fortran and Python. Python loads an (n1,n2,n3) Fortran matrix as: (n3,n2,n1),
# and this script keeps that consistent so that the final H5 datasets are really
# what would be saved by CUDA Fortran. This allows QCR visualization scripts to
# be used on CUDA QCR data directly.

import h5py
import numpy as np
from scipy.io import FortranFile

# set CUDA QCR version
cver = '2d1e'

# specify precision of saved data
dtype = 'float32'

# n = (nx, ny, nt, ntg)
n = np.loadtxt('dims.dat').astype(int)

# time
f = FortranFile('t.dat','r')
t = np.zeros((1,n[2]), dtype=dtype)
t[0,:] = np.array(f.read_reals(dtype=dtype))
t = np.transpose(t)
# time, grid
f = FortranFile('tg.dat','r')
tg = np.zeros((1,n[3]), dtype=dtype)
tg[0,:] = np.array(f.read_reals(dtype=dtype))
tg = np.transpose(tg)
# x/x1 position
if (cver == '2d1e'): f = FortranFile('x.dat','r')
if (cver == '1d2e'): f = FortranFile('x1.dat','r')
x = np.zeros((1,n[0]), dtype=dtype)
x[0,:] = np.array(f.read_reals(dtype=dtype))
x = np.transpose(x)
# y/x2 position
if (cver == '2d1e'): f = FortranFile('y.dat','r')
if (cver == '1d2e'): f = FortranFile('x2.dat','r')
y = np.zeros((1,n[1]), dtype=dtype)
y[0,:] = np.array(f.read_reals(dtype=dtype))
y = np.transpose(y)
# absorbing potential
f = FortranFile('Va.dat','r')
Va = np.array(f.read_reals(dtype=dtype)).reshape(n[1],n[0])
# electric field
if (cver == '2d1e'):
    Ert = np.zeros((2,n[2]), dtype=dtype)
    f = FortranFile('Ex.dat','r')
    Ert[0,:] = np.array(f.read_reals(dtype=dtype))
    f = FortranFile('Ey.dat','r')
    Ert[1,:] = np.array(f.read_reals(dtype=dtype))
    Ert = np.transpose(Ert)
elif (cver == '1d2e'):
    Ert = np.zeros((1,n[2]), dtype=dtype)
    f = FortranFile('Et.dat','r')
    Ert[0,:] = np.array(f.read_reals(dtype=dtype))
    Ert = np.transpose(Ert)
# vector potential
if (cver == '2d1e'):
    Art = np.zeros((2,n[2]), dtype=dtype)
    f = FortranFile('Ax.dat','r')
    Art[0,:] = np.array(f.read_reals(dtype=dtype))
    f = FortranFile('Ay.dat','r')
    Art[1,:] = np.array(f.read_reals(dtype=dtype))
    Art = np.transpose(Art)
elif (cver == '1d2e'):
    Art = np.zeros((1,n[2]), dtype=dtype)
    f = FortranFile('At.dat','r')
    Art[0,:] = np.array(f.read_reals(dtype=dtype))
    Art = np.transpose(Art)
# wavefunction
f = FortranFile('Re_psi.dat','r')
Re_psi = np.array(f.read_reals(dtype=dtype)).reshape(n[3],n[1],n[0])
f = FortranFile('Im_psi.dat','r')
Im_psi = np.array(f.read_reals(dtype=dtype)).reshape(n[3],n[1],n[0])
# time-dependent potential
f = FortranFile('Vp.dat','r')
Vp = np.array(f.read_reals(dtype=dtype)).reshape(n[3],n[1],n[0])
# norm psi
f = FortranFile('norm_psi.dat','r')
norm_psi = np.zeros((1,n[2]), dtype=dtype)
norm_psi[0,:] = np.array(f.read_reals(dtype=dtype))
norm_psi = np.transpose(norm_psi)
# energy psi
f = FortranFile('energy_psi.dat','r')
energy_psi = np.zeros((1,n[2]), dtype=dtype)
energy_psi[0,:] = np.array(f.read_reals(dtype=dtype))
energy_psi = np.transpose(energy_psi)

hf = h5py.File('qm_data.h5','w')
hf.create_dataset('t', data=t)
hf.create_dataset('tg', data=tg)
if (cver == '2d1e'):
    hf.create_dataset('x', data=x)
    hf.create_dataset('y', data=y)
elif (cver == '1d2e'):
    hf.create_dataset('x1', data=x)
    hf.create_dataset('x2', data=y)
hf.create_dataset('Va', data=Va)
if (cver == '2d1e'):
    hf.create_dataset('Ert', data=Ert)
    hf.create_dataset('Art', data=Art)
elif (cver == '1d2e'):
    hf.create_dataset('Et', data=Ert)
    hf.create_dataset('At', data=Art)
hf.create_dataset('Re(psi)', data=Re_psi)
hf.create_dataset('Im(psi)', data=Im_psi)
hf.create_dataset('Vp', data=Vp)
hf.create_dataset('norm_psi', data=norm_psi)
hf.create_dataset('energy_psi', data=energy_psi)
hf.close()
