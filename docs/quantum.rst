Documentation
=============
Module: Quantum Mechanics ``(quantum.f08)``
-------------------------------------------
Revised 9/14/2025

**Classes**

.. code-block::
    :number-lines:

    type SchrodingerWavefunction2D
      Schrodinger 2D Wavefunction class
    variables
      psi: wavefunction data
      grad_psi: wavefunction gradient
      phase: wavefunction phase
      norm: normalization value
      energy: Hamiltonian expectation value
    variables for propagation
      D2(x,y), M2(x,y): sparse Crank-Nicolson-Numerov (CNN) 2nd derivative matrices
      stab(x,y): row/column population table for sparse tri-diagonal matrices
    procedures
      init_vars => psi2d_initialize_vars: initialize the wavefunction object and its data arrays
      init_form => psi2d_initialize_form: generate a random wavefunction of definite parity
      make_cnn_mats => psi2d_make_cnn_matrices: initialize the sparse Crank-Nicolson-Numerov derivative matrices
      propagate_fft => psi2d_dt_propagate_fft: advance the wavefunction by dt using the split-operator method
      propagate_cnn => psi2d_dt_propagate_cnn: advance the wavefunction by dt using the Crank-Nicolson-Numerov method
      destroy => psi2d_destructor: object destructor
    note
      - The variable dimensions are:
        psi, phase: (nx,ny)
        grad_psi: (nx,ny,2)
        norm, energy: (1,nt)
        D2(x,y), M2(x,y): (1,3*n(x,y)-2)
        stab(x,y): (3*n(x,y)-2,2)

.. code-block::
    :number-lines:

    type SchrodingerWavefunction1D
      Schrodinger 1D Wavefunction class
    variables, variables for propagation
      See 2D Wavefunction class
    procedures
      init_vars => psi1d_initialize_vars: initialize the wavefunction object and its data arrays
      init_form => psi1d_initialize_form: generate a random wavefunction of definite parity
      make_cnn_mats => psi1d_make_cnn_matrices: initialize the sparse Crank-Nicolson-Numerov derivative matrices
      propagate_fft => psi1d_dt_propagate_fft: advance the wavefunction by dt using the split-operator method
      propagate_cnn => psi1d_dt_propagate_cnn: advance the wavefunction by dt using the Crank-Nicolson-Numerov method
      destroy => psi1d_destructor: object destructor
    note
      - The variable dimensions are:
        psi, grad_psi, phase: (1,nx)
        norm, energy: (1,nt)
        D2x, M2x: (1,3*nx-2)
        stab: (3*nx-2,2)

.. code-block::
    :number-lines:

    type SchrodingerWavefunction1DR
      Schrodinger 1D Radial Wavefunction class
    variables
      phi: matrix of radial functions
      V0: atomic potential
      dV0_dr: radial derivative of V0
      norm: normalization value
      energy: Hamiltonian expectation value
      Z, m: nuclear charge, magnetic quantum number
    variables for propagation
      Va: absorbing potential
      D1r, D2r: sparse Crank-Nicolson derivative matrices
      M1r, M2r: sparse Muller matrices
      D211, M211: upper-element matrix corrections for l=m=0
      stab: row/column population table for sparse tri-diagonal matrices
      clm: orbital angular momentum matrix elements
    procedures
      init_vars => psi1dr_initialize_vars: initialize the wavefunction object and its data arrays
      init_prop => psi1dr_initialize_propagators: initialize all propagation matrices
      init_form => psi1dr_initialize_form: generate a random initial wavefunction
      propagate => psi1dr_dt_propagate_free, psi1dr_dt_propagate_full: advance the wavefunction by dt
      prep_atom => psi1dr_prepare_atomic_state: prepare the wavefunction in a bound atomic eigenstate
      destroy => psi1dr_destructor: object destructor
    note
      - The variable dimensions are:
        phi, V0, dV0_dr: (l_max+1,nr)
        norm, energy: (1,nt)
        D1/2r, M1/2r: (1,3*nr-2)
        stab: (3*nr-2,2)
        clm: (1,l_max)
        Va: (1,nr)
      - The magnetic quantum number (m) is fixed by the initial state.
      - The absorbing potential (Va) must be initialized externally.
      - The prep_atom subroutine works via Imaginary Time Propagation (ITP), it can only produce states
        for which n=l+1, where (n,l) are the principal and angular momentum quantum numbers, respectively.
        For the general case, use the Eigenstate Distillation Method (EDM) as described in:
        D. Bauer, Computational Strong-Field Quantum Dynamics, Chap. II, Sec. 2.2.5.

.. code-block::
    :number-lines:

    type tSURFF2D
      Time-dependent surface-flux method (2D) class
    variables
      p_dist: probability amplitude of the photo-electron momentum distribution
      kx, ky: momentum distribution bins
      xl, yl: discrete surface integration points
      ixl, iyl: lower-bound coordinate grid indices
      interp: wavefunction interpolation method
        (nn: nearest-neighbor, b3: bicubic)
      R0, dphi: surface radius, azimuthal angle step-size
      dti, iti: integration time period, time-step period
      Ns: no. surface points
      nk: no. distribution bins
      k(x,y)_lim: distribution extents
      enable: toggle tSURFF calculation
    procedures
      init => tSURFF2D_initialize: initialize the tSURFF2D object and its data arrays
      dt_step => tSURFF2D_dt_step: dti-surface-integrate the probability amplitude
      destroy => tSURFF2D_destructor: object destructor

.. code-block::
    :number-lines:

    type pconst_mks
      Fundamental physical constants, meter-kilogram-second (MKS) base units.
      NIST CODATA 2018 recommended values.

.. code-block::
    :number-lines:

    type pconst_cgs
      Fundamental physical constants, centimeter-gram-second (CGS) base units.
      NIST CODATA 2018 recommended values.

**Subroutines** ``(initialize_vars)``

.. code-block::
    :number-lines:

    subroutine psi1d_initialize_vars
      Initialize a 1D wavefunction object.
    input
      this: Schrodinger 1D wavefunction class object
      nx, nt: no. spatial/temporal grid points

.. code-block::
    :number-lines:

    subroutine psi1dr_initialize_vars
      Initialize a radial 1D wavefunction object.
    input
      this: Schrodinger 1DR wavefunction class object
      nr, nt: no. spatial/temporal grid points
      l_max: azimuthal quantum numbers in expansion

.. code-block::
    :number-lines:

    subroutine psi2d_initialize_vars
      Initialize a 2D wavefunction object.
    input
      this: Schrodinger 2D wavefunction class object
      nr, nt: no. spatial/temporal grid points

**Subroutines** ``(destructors)``

.. code-block::
    :number-lines:

    subroutine psi1d_destructor, psi1dr_destructor, psi2d_destructor
      Deallocates the wavefunction object.
    input
      this: Schrodinger 1D, 1DR, or 2D wavefunction class object

.. code-block::
    :number-lines:

    subroutine tSURFF2D_destructor
      Deallocates the tSURFF2D object.
    input
      this: tSURFF2D class object

**Subroutines** ``(initialize_form)``

.. code-block::
    :number-lines:

    subroutine psi1d_initialize_form
      Create a random 1D wavefunction of definite parity.
    input
      this: Schrodinger 1D wavefunction class object
      parity: desired eigenstate parity (+/-1)
      dx: spatial step-size

.. code-block::
    :number-lines:

    subroutine psi1dr_initialize_form
      Create a random 1D radial wavefunction.
    input
      this: Schrodinger 1DR wavefunction class object
      dr: spatial step-size

.. code-block::
    :number-lines:

    subroutine psi2d_initialize_form
      Create a random 2D wavefunction of definite parity.
    input
      this: Schrodinger 2D wavefunction class object
      parity: desired eigenstate parity (+/-1)
      dr: spatial step-size

**Subroutines** ``(make_cnn_matrices)``

.. code-block::
    :number-lines:

    subroutine psi1d_make_cnn_matrices
      Create the sparse Crank-Nicolson-Numerov 2nd derivative matrices.
    input
      this: Schrodinger 1D wavefunction class object
      nx, dx: no. spatial grid points/step-size
    output
      D2x, M2x: sparse 2nd derivative matrices
      stab: table of sparse indices

.. code-block::
    :number-lines:

    subroutine psi2d_make_cnn_matrices
      Create the sparse Crank-Nicolson-Numerov 2nd derivative matrices.
    input
      this: Schrodinger 2D wavefunction class object
      nr, dr: no. spatial grid points/step-size
    output
      D2(x,y), M2(x,y): sparse 2nd derivative matrices
      stab(x,y): table of sparse indices

**Subroutines** ``(initialize_propagators)``

.. code-block::
    :number-lines:

    subroutine psi1dr_initialize_propagators
      Create all propagation matrices for a 1D radial wavefunction.
    input
      this: Schrodinger 1DR wavefunction class object
      r, dr: spatial grid/step-size
    output
      V0, dV0_dr: atomic potential with centrifugal term
      D1/2r, M1/2r: r-space derivative matrices
      D211, M211: D/M2r upper-element corrections
      stab: table of sparse indices
      clm: angle-space rotation matrix elements

**Subroutines** ``(dt_propagate)``

.. code-block::
    :number-lines:

    subroutine psi1d_dt_propagate_fft
      Advance the Schrodinger wavefunction (psi) by dt using the split-operator method.
    input
      this: Schrodinger 1D wavefunction class object
      T, V: kinetic and potential energy arrays of dim(1,nx)
      dx, dt: spatial/temporal step-size
      j: imaginary time constant
    output
      psi: updated wavefunction; psi(t+dt)
    note
      - Assumes V = V(x), time-independent.
      - The potential is generally complex.
      - The wavefunction is over-written.
      - Pass j = -i for Imaginary Time Propagation (ITP) and cmplx(1.0) otherwise.
      - ITP results in non-unitary dynamics, so normalization must be enforced manually.

.. code-block::
    :number-lines:

    subroutine psi1d_dt_propagate_cnn
      Advance the Schrodinger wavefunction (psi) by dt using the Crank-Nicolson-Numerov method.
    input
      this: Schrodinger 1D wavefunction class object
      V: potential energy array of dim(1,nx)
      dt: temporal step-size
      j: imaginary time constant
    output
      psi: updated wavefunction; psi(t+dt)
    note
      - Assumes V = V(x), time-independent.
      - The potential is generally complex.
      - The wavefunction is over-written.
      - Pass j = -i for Imaginary Time Propagation (ITP) and cmplx(1.0) otherwise.
      - ITP results in non-unitary dynamics, so normalization must be enforced manually.

.. code-block::
    :number-lines:

    subroutine psi1dr_prepare_atomic_state
      Propagate one radial eigenfunction of the Schrodinger state (phi) through imaginary time.
      Used to prepare an initial atomic (l,m)-state.
    input
      this: Schrodinger 1DR wavefunction class object
      lp: desired OAM state to distill
      ntau: no. relaxation time-steps
      dr, dt: spatial/temporal step-size
      pure: logical, kill every (l,m) component except lp
    output
      phi: distilled atomic wavefunction
    note
      - The wavefunction is over-written.
      - ITP renormalization is performed in this routine.

.. code-block::
    :number-lines:

    subroutine psi1dr_dt_propagate_free
      Advance the radial eigenfunctions of the Schrodinger state (phi) by dt.
      Atomic potential only; No external electromagnetic fields.
    input
      this: Schrodinger 1DR wavefunction class object
      dt: temporal step-size
    output
      phi: updated wavefunction; phi(t+dt)
    note
      - The wavefunction is over-written.
      - The complex absorbing potential is used.

.. code-block::
    :number-lines:

    subroutine psi1dr_dt_propagate_full
      Advance the radial eigenfunctions of the Schrodinger state (phi) by dt.
      Time-propagation with a linearly-polarized field.
    input
      this: Schrodinger 1DR wavefunction class object
      A: instantaneous field vector potential
      r: spatial grid
      dt: temporal step-size
    output
      phi: updated wavefunction; phi(t+dt)
    note
      - The wavefunction is over-written.
      - The intermediate r-space transformation calls on the _free propagation routine.

.. code-block::
    :number-lines:

    subroutine psi2d_dt_propagate_fft
      Advance the Schrodinger wavefunction (psi) by dt using the split-operator method.
    input
      this: Schrodinger 2D wavefunction class object
      T, V: kinetic and potential energy arrays of dim(nr(1),nr(2))
      dr, dt: spatial/temporal step-size
      j: imaginary time constant
    output
      psi: updated wavefunction; psi(t+dt)
    note
      - Assumes V = V(x,y), time-independent.
      - The potential is generally complex.
      - The wavefunction is over-written.
      - Pass j = -i for Imaginary Time Propagation (ITP) and cmplx(1.0) otherwise.
      - ITP results in non-unitary dynamics, so normalization must be enforced manually.

.. code-block::
    :number-lines:

    subroutine psi2d_dt_propagate_cnn
      Advance the Schrodinger wavefunction (psi) by dt using the Crank-Nicolson-Numerov method.
    input
      this: Schrodinger 2D wavefunction class object
      V: potential energy array of dim(nr(1),nr(2))
      dt: temporal step-size
      j: imaginary time constant
    output
      psi: updated wavefunction; psi(t+dt)
    note
      - Assumes V = V(x,y), time-independent.
      - The potential is generally complex.
      - The wavefunction is over-written.
      - Pass j = -i for Imaginary Time Propagation (ITP) and cmplx(1.0) otherwise.
      - ITP results in non-unitary dynamics, so normalization must be enforced manually.

**Functions** ``(expectE)``

.. code-block::
    :number-lines:

    function expectE_ND_fft
      Calculates the energy expectation value given the ND wavefunction.
      FFT version.
    input
      this: Schrodinger ND wavefunction class object
      dr, dp: space and momentum step sizes
      T, V: kinetic and potential energy arrays of dim(nr)
    note
      - The input potential must be real.
      - Every dimension of psi must be a power of 2.

.. code-block::
    :number-lines:

    function expectE_ND_cnn
      Calculates the energy expectation value given the ND wavefunction.
      Crank-Nicolson version.
    input
      this: Schrodinger ND wavefunction class object
      V: potential energy array of dim(nr)
      dr: spatial step-size
    note
      - The input potential must be real.

.. code-block::
    :number-lines:

    function expectE_1DR
      Calculates the energy expectation value given the 1D radial wavefunction.
    input
      this: Schrodinger 1DR wavefunction class object
      dr: spatial step-size

**Subroutines** ``(photoe_spectrum_winop)``

.. code-block::
    :number-lines:

    subroutine photoe_spectrum_winop_1D
      Calculates the photo-electron spectrum using
      the nth-order energy window operator method.
      1D Cartesian wavefunction.
    input
      psi: 1D wavefunction
      V0: atomic potential
      E: energy bins
      dx: spatial step-size
      n: window order
    output
      W: spectrum

.. code-block::
    :number-lines:

    subroutine photoe_spectrum_winop_1DR
      Calculates the photo-electron spectrum using
      the nth-order energy window operator method.
      1D Radial wavefunction.
    input
      this: Schrodinger 1DR wavefunction class object
      E: energy bins
      dr: spatial step-size
      n: window order
    output
      W: spectrum

**Subroutines** ``(tSURFF2D)``

.. code-block::
    :number-lines:

    subroutine tSURFF2D_initialize
      Initialize the tSURFF2D object and its data arrays.
      The user must externally supply values for the following class parameters:
        enable, interp, dti, Ns, R0, k(x,y)_lim, nk
    input
      S: tSURFF2D class object
      x, y: spatial mesh arrays
      dt: global time-step
    note
      - The spatial indices (ixl,iyl) correspond to the nearest upper-left point
        in the underlying spatial mesh (spanned by inputs x & y).
      - The integration time-step period (iti) is an integer, given by the
        nearest whole ratio of dti to the global time-step.

.. code-block::
    :number-lines:

    subroutine tSURFF2D_dt_step
      Perform a surface-flux integration, advancing in time by dti.
    input
      S: tSURFF2D class object
      wavefn: Schrodinger 2D wavefunction class object
      A, C: two-component vector potential and excursion at the current time (t)
      x, y: spatial mesh arrays
      t: the current time
    output
      p_dist: dti-advanced momentum probability amplitude
    note
      - A simpler 1D version of the underlying theory is developed in:
          D. Bauer, Computational Strong-Field Quantum Dynamics.
        See also:
          V. Mosert and D. Bauer, "Photoelectron spectra with Qprop and t-SURFF",
            Comput. Phys. Commun. 207, 452 (2016).
      - The wavefunction at the current time-step is interpolated via either:
          1. nearest-neighbor (nn), approximating psi(R0,t) by the value at the nearest upper-left grid point.
        or
          2. bicubic (b3), interpolating the 4 surrounding mesh points to compute psi(R0,t).
        Generally, b3 is far more computationally demanding.

**Subroutines** ``(radiative_intensity)``

.. code-block::
    :number-lines:

    note
      - The radiation spectrum is the Fourier transform of the intensity time-series array S(t),
        with the frequency calculated from the temporal grid: ω = 2π*fft_freq(nt,dt).
      - See D. Bauer, Computational Strong-Field Quantum Dynamics, Chap. II, Sec. 3.

.. code-block::
    :number-lines:

    subroutine radiative_intensity_1DR
      Calculates the electron radiation emission intensity at the current time.
      1D Radial wavefunction.
    input
      this: Schrodinger 1DR wavefunction class object
      Et: instantaneous electric field strength
      dr: spatial step-size
    output
      S: instantaneous radiative intensity

.. code-block::
    :number-lines:

    subroutine radiative_intensity_2D1e
      Calculates the electron radiation emission intensity at the current time.
      2D wavefunction (1 two-dimensional electron).
    input
      this: Schrodinger 2D wavefunction class object
      aV: matrix of partial derivatives of the bare atomic potential, dim(nx,ny,2)
      Et: instantaneous electric field strength (two polarization components)
      dr: spatial step-size
    output
      S: instantaneous radiative intensity

.. code-block::
    :number-lines:

    subroutine radiative_intensity_1D2e
      Calculates the electron radiation emission intensity at the current time.
      2D wavefunction (2 one-dimensional electrons).
    input
      this: Schrodinger 2D wavefunction class object
      aV: matrix of partial derivatives of the bare atomic potential, dim(nx,ny,2)
      Et: instantaneous electric field strength
    output
      S: instantaneous radiative intensity
    note
      - S(1,2) corresponds to electrons 1 & 2, respectively.
      - The total emission spectrum can be obtained by Fourier-transforming S(1)(t) + S(2)(t).

**Subroutines** ``(bohm)``

.. code-block::
    :number-lines:

    subroutine calc_bohm_velocity
      Compute the Bohmian velocity at the current time using the probability current.
    input
      nx, nt, k: no. spatial/temporal grid points and current time-step index
      dx: spatial step-size
      psi: wavefunction, dim(nt,nx)
    output
      bv: Bohmian velocity; v(x,t) = J(x,t)/rho(x,t) where rho(x,t) = abs(psi(x,t))^2
      Jxt: probability current; J(x,t)
    note
      - Outputs have dim(nt,nx).

.. code-block::
    :number-lines:

    subroutine calc_bohm_velocity_from_phase
      Compute the Bohmian velocity at the current time using the phase of the wavefunction, S(x,t).
    input
      nx, nt, k: no. spatial/temporal grid points and current time-step index
      dx: spatial step-size
      psi: wavefunction, dim(nt,nx)
    output
      bv: Bohmian velocity; v(x,t) = grad(S(x,t))
    note
      - This routine is NOT recommended over calc_bohm_velocity due to the phase-unwrapping problem.

.. code-block::
    :number-lines:

    subroutine calc_bohm_trajectories
      Obtain Bohmian trajectories from pre-computed velocity field information.
    input
      nx, nt: no. spatial/temporal grid points
      dt: temporal step-size
      x: spatial grid
      bv: Bohm velocity array, dim(nt,nx)
    output
      bx: Bohmian trajectory array
    note
      - Each trajectory is initialized by the first row of bv, i.e., v(x,t=0).

**Subroutines** ``(misc)``

.. code-block::
    :number-lines:

    subroutine chk_continuity_eqn
      Evaluate the continuity equation for all time.
    input
      nx, nt: no. spatial/temporal grid points and current time-step index
      dx, dt: spatial/temporal step-size
      psi, Jxt: wavefunction and probability current, dim(nt,nx)
    output
      cty: evaluated continuity equation

**Functions** ``(misc)``

.. code-block::
    :number-lines:

    function E_hydrogen
      Exact hydrogen energy levels in atomic units.
    input
      n, l: principal, azimuthal quantum number
    output
      En: energy
    note
      - See S. Weinberg, Quantum Theory of Fields, Vol. I, Eq. (1.1.27).
