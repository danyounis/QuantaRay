![](docs/logos/v2-small.png)

A library for strong-field quantum dynamical calculations.

### Programs
**`main`**: Extended virtual detector method — the time-dependent Schrödinger equation \(TDSE\) is integrated on a fixed lattice close to an atomic nucleus, and information is encoded within Bohmian-like particles which are subsequently propagated beyond the absorbing boundary \[1–4\].

**`prstrm`**: PairStorm — A quantum electrodynamics (QED) Monte Carlo particle code. Stochastic model of high-energy photon emission and electron-positron pair production in a prescribed background field, within the framework of the Volkov dressed-state representation and local-constant-field approximation \[5\].

**`radse`**: A pseudo-spectral representation of the radial Schrödinger wavefunction, based on the Patchkovskii–Muller propagator \[6\].

**`rabi-split`**: Similar to `radse`, this code is used to model atomic energy-level splitting in an electromagnetic field. The Kulander–Schafer window operator \[7\] is employed to calculate the photoelectron spectrum \[8\].

**`cuda`**: A test implementation of the NVIDIA CUDA® library to solve the TDSE. *\(Under development\)*

### Modules
**`cufft.cuf`**: NVIDIA CUDA® Fast Fourier Transform (FFT) interface

**`emfm.f08`**: Electromagnetic field module

**`math.f08`**: Core mathematical functions

**`optimize.f08`**: Numerical optimization procedures

**`prec.F08`**: Numerical precision types

**`qed_mcpm.F08`**: QED Monte Carlo particle module

**`quantum.f08`**: Quantum mechanics module

**`rochester.f08`**: Rochester potential functions

**`vdm.f08`**: Virtual detector module

### Dependencies
BlueHive HPC modules:
 - `gcc/11.2.0/b1`
 - `nvhpc/20.7`
 - `openmpi/2.1.1/b1`
 - `hdf5/1.12.1/b1`
 - `anaconda3/2021.11`
 - `lapack/3.9.0/b2`
 - `openblas/0.3.10/b1`

MPI-specific:
 - `gcc/9.1.0`
 - `openmpi/4.0.4/b4`

### References
\[1\] B. Feuerstein and U. Thumm, *"On the computation of momentum distributions within wavepacket propagation calculations"*, [J. Phys. B **36**, 707 \(2003\)](https://doi.org/10.1088/0953-4075/36/4/305).

\[2\] X. Wang, J. Tian, and J.H. Eberly, *"Extended Virtual Detector Theory for Strong-Field Atomic Ionization"*, [Phys. Rev. Lett. **110**, 243001 \(2013\)](https://doi.org/10.1103/PhysRevLett.110.243001).

\[3\] R.-H. Xu and X. Wang, *"Extended virtual detector theory including quantum interferences"*, [AIP Adv. **11**, 025124 \(2021\)](https://doi.org/10.1063/5.0040193).

\[4\] D. Younis and J.H. Eberly, *"Strong-field nonsequential double photoionization using virtual-detector theory with path summation"*, [Phys. Rev. A **107**, 053117 \(2023\)](https://doi.org/10.1103/PhysRevA.107.053117).

\[5\] C.P. Ridgers *et al.*, *"Modelling gamma-ray photon emission and pair production in high-intensity laser–matter interactions"*, [J. Comp. Phys. **260**, 273 \(2014\)](https://doi.org/10.1016/j.jcp.2013.12.007).

\[6\] S. Patchkovskii and H.G. Muller, *"Simple, accurate, and efficient implementation of 1-electron atomic time-dependent Schrödinger equation in spherical coordinates"*, [Comput. Phys. Commun. **199**, 153 \(2016\)](https://doi.org/10.1016/j.cpc.2015.10.014).

\[7\] K.J. Schafer and K.C. Kulander, *"Energy analysis of time-dependent wave functions: Application to above-threshold ionization"*, [Phys. Rev. A **42**, 5794\(R\) \(1990\)](https://doi.org/10.1103/PhysRevA.42.5794).

\[8\] D. Younis and J.H. Eberly, *"Benchmark of few-level quantum theory... for the strong-field Autler–Townes"*, [J. Phys. B **55**, 164001 \(2022\)](https://doi.org/10.1088/1361-6455/ac7d7f).
