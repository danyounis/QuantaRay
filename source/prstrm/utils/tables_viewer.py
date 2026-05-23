# tables_viewer.py
import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
font={'size':10,'family':'monospace'}; lw=1; rc('font',**font);

'''-----------------------------------------------------------------------
Synchrotron function, to determine the γ emission rate as a function of η.
-----------------------------------------------------------------------'''

with open('../tables/eta_2i.bin','rb') as fi:
    _ = np.fromfile(fi, dtype=np.int32, count=1)
    eta_2i = np.fromfile(fi, dtype=float)

with open('../tables/h_2i.bin','rb') as fi:
    _ = np.fromfile(fi, dtype=np.int32, count=1)
    h_2i = np.fromfile(fi, dtype=float)

plt.semilogx(eta_2i, h_2i, lw=lw, c='k')
ax = plt.gca()
ax.set_xlabel('Lepton quantum parameter, η')
ax.set_ylabel('Synchrotron function, h(η)')
ax.set_title('Files: eta_2i.bin, h_2i.bin', fontsize=font['size'])
plt.show()

'''----------------------------------------------------------------------------
Emissivity function, to determine the e-/e+ production rate as a function of χ.
----------------------------------------------------------------------------'''

with open('../tables/chi_3i.bin','rb') as fi:
    _ = np.fromfile(fi, dtype=np.int32, count=1)
    chi_3i = np.fromfile(fi, dtype=float)

with open('../tables/g_3i.bin','rb') as fi:
    _ = np.fromfile(fi, dtype=np.int32, count=1)
    g_3i = np.fromfile(fi, dtype=float)

plt.semilogx(chi_3i, g_3i, lw=lw, c='k')
ax = plt.gca()
ax.set_xlabel('Photon quantum parameter, χ')
ax.set_ylabel('Emissivity function, g(χ)')
ax.set_title('Files: chi_3i.bin, g_3i.bin', fontsize=font['size'])
plt.show()

'''--------------------------------------------------------------------
Tables to determine the probability of an η-lepton emitting a χ-photon.
--------------------------------------------------------------------'''

with open('../tables/eta_0i.bin','rb') as fi:
    _ = np.fromfile(fi, dtype=np.int32, count=1)
    eta_0i = np.fromfile(fi, dtype=float)

with open('../tables/chi_0ji.bin','rb') as fi:
    shape = np.fromfile(fi, dtype=np.int32, count=2)
    chi_0ji = np.fromfile(fi, dtype=float).reshape(shape, order='F')

with open('../tables/Py_0ij.bin','rb') as fi:
    shape = np.fromfile(fi, dtype=np.int32, count=2)
    Py_0ij = np.fromfile(fi, dtype=float).reshape(shape, order='F')

fig, ax = plt.subplots(1,3, figsize=(16,5), dpi=100, constrained_layout=True)

ax[0].semilogy(eta_0i, lw=lw, c='k')
ax[0].set_xlabel('Row index i')
ax[0].set_ylabel('Lepton quantum parameter, η(i)')
ax[0].set_title('File: eta_0i.bin', fontsize=font['size'])

im = ax[1].imshow(np.log10(chi_0ji), aspect=0.5, origin='lower', clim=[-6,4], cmap='binary')
ax[1].set_xlabel('Column index j')
ax[1].set_ylabel('Row index i')
ax[1].set_title('File: chi_0ji.bin', fontsize=font['size'])
cbar = plt.colorbar(im, shrink=0.95)
cbar.set_label('Element (i,j) = log χ(i,j) corresponding to η(i)', fontsize=font['size'])

im = ax[2].imshow(Py_0ij, aspect=0.5, origin='lower', clim=[0,1], cmap='RdBu_r')
ax[2].set_xlabel('Column index j')
ax[2].set_ylabel('Row index i')
ax[2].set_title('File: Py_0ij.bin', fontsize=font['size'])
cbar = plt.colorbar(im, shrink=0.95)
cbar.set_label('''Element (i,j) = Probability of
emitting χ(i,j) photon by η(i) lepton''', fontsize=font['size'])

plt.show()

'''--------------------------------------------------------------------
Tables to determine the energy split of a χ-photon between e-/e+ pairs.
--------------------------------------------------------------------'''

with open('../tables/chi_1i.bin','rb') as fi:
    _ = np.fromfile(fi, dtype=np.int32, count=1)
    chi_1i = np.fromfile(fi, dtype=float)

with open('../tables/frac_1j.bin','rb') as fi:
    _ = np.fromfile(fi, dtype=np.int32, count=1)
    frac_1j = np.fromfile(fi, dtype=float)

with open('../tables/Pf_1ij.bin','rb') as fi:
    shape = np.fromfile(fi, dtype=np.int32, count=2)
    Pf_1ij = np.fromfile(fi, dtype=float).reshape(shape, order='F')

fig, ax = plt.subplots(1,3, figsize=(16,4.5), dpi=100, constrained_layout=True)

ax[0].semilogy(chi_1i, lw=lw, c='k')
ax[0].set_xlabel('Row index i')
ax[0].set_ylabel('Photon quantum parameter, χ(i)')
ax[0].set_title('File: chi_1i.bin', fontsize=font['size'])

ax[1].plot(frac_1j, lw=lw, c='k')
ax[1].set_xlabel('Column index j')
ax[1].set_ylabel('Energy-share fraction, f(j)')
ax[1].set_title('File: frac_1j.bin', fontsize=font['size'])

im = ax[2].imshow(Pf_1ij, origin='lower', clim=[0,1], cmap='RdBu_r')
ax[2].set_xlabel('Column index j')
ax[2].set_ylabel('Row index i')
ax[2].set_title('File: Pf_1ij.bin', fontsize=font['size'])
cbar = plt.colorbar(im)
cbar.set_label('''Element (i,j) = Probability of a χ(i) photon giving
fractions f(j) & 1-f(j) of its energy to the e- & e+''', fontsize=font['size']-1)

plt.show()
