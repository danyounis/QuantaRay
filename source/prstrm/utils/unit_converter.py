# unit_converter.py
import numpy as np

class cgs:
    c = 2.99792458e+10 # speed of light [cm/s] (exact)
    e = 4.803204712570263e-10 # elementary charge [statC] (exact)
    m_e = 9.1093837015e-28 # electron mass [g] (+/- 2.8e-37 g)
    hbar = 1.0545718176461565e-27 # reduced Planck constant [erg.s] (exact)
    alpha = 7.2973525693e-3 # fine-structure constant (+/- 1.1e-12)
    lambda_Compton = 2.42631023867e-10 # Compton wavelength [cm] (+/- 7.3e-20 cm)

'''
Units, multiply for [sim.u.]->[CGS]
    [length] = 1/k0_icm, [mass] = m_e, [time] = 1/w0_Hz, [charge] = |e|,
    [speed] = w0_Hz/k0_icm = c, [momentum] = m_e*c, [energy] = m_e*(c^2),
    [field] = w0_Hz*m_e*c/|e|

Notes
- In expressions, m_e, c, and |e| are unity, while ℏ becomes the `hbar_` variable.
  For example, the QED critical field Ec [CGS] = (m_e^2)*(c^3)/(|e|*ℏ) is simply 1/hbar_ in [sim.u.].
- hbar_ and q_ are related by the fine-structure constant.
'''
(k0_, w0_, lam0_) = (1.0, 1.0, 2*np.pi)
# do not modify the above line; it is for reference

'''-------
Parameters
-------'''
lam0_cm = 1.e-4 # normalization length [cm]

'''---------------
Derived Quantities
---------------'''
k0_icm = (2*np.pi)/lam0_cm # wavenumber [1/cm]
w0_Hz = cgs.c*k0_icm # frequency [Hz]
hbar_ = cgs.hbar*k0_icm/cgs.m_e/cgs.c # reduced Planck constant [sim.u.]
q_ = np.sqrt(cgs.alpha*hbar_) # electric charge [sim.u.]
F_ = (w0_Hz*cgs.m_e*cgs.c)/cgs.e # E/B-field conversion factor [CGS]->[sim.u.]
