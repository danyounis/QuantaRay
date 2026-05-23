Documentation
=============
Module: Virtual Detector ``(vdm.f08)``
--------------------------------------
Revised 9/14/2025

**Classes**

.. code-block::
    :number-lines:

    type vdet
      Virtual Detector (VD) class
    parameters
      xl, yl: VD position coordinates
      ixl, iyl: lower-bound coordinate grid indices
    variables
      Krt: recorded momentum
      Jrt: recorded probability current density
      rho: recorded probability density
      phase: recorded wavefunction phase
    procedures
      init => vdet_initialize: initialize a virtual detector and its data arrays.
      trigger(1/2) => vdet_calculate_current_(1/2)d: calculate and record the instantaneous momentum, probability current, probability density, and phase.
    note
      - The variables Krt and Jrt have dimensions (nr,nt), where nr is the number of
        spatial components and nt is the number of time steps; and phase and rho have dimensions (1,nt).

.. code-block::
    :number-lines:

    type edet
      End Detector (ED) class
    variables
      nde: no. detected electrons
      bfwt: bound/free weight totals
      data: recorded electron information
    procedures
      escan: record electron information
    note
      - The data variable has dimensions (nde,4) (1D) or (nde,6) (2D).
      - An ED entry contains electron trajectory:
        (x, px, phase, weight) (1D) or (x, y, px, py, phase, weight) (2D).
      - bfwt(1,2) = sum of (bound,free) virtual electron weights

.. code-block::
    :number-lines:

    type particle_electron
      Electron class
    variables
      x, y: position coordinates
      px, py: momentum components
      weight: statistical weight
      phase: accumulated trajectory phase
      ix, iy: position grid indices
      propagate: toggle electron dynamics
    procedures
      apush(1/2) => electron_dt_propagate_analytic_(1/2)d: analytic propagation routine
      npush => electron_dt_propagate_numeric: numeric propagation routine

.. code-block::
    :number-lines:

    type grid_electron
      Electron Grid class
    parameters
      x_lim, y_lim: domain boundaries
      dr, nr: step size and no. points
    variables
      x, y: spatial mesh arrays

**Subroutines** ``(initialize)``

.. code-block::
    :number-lines:

    subroutine vdet_initialize
      Initialize a virtual detector.
    input
      this: virtual detector class object
      geom: detector geometry in x-y space
      R0: detector shape characteristic size
      Nv: total no. virtual detectors
      x, y: spatial mesh arrays
      n: this virtual detector number
      nt: no. temporal grid points
      sdim: dimensionality of the simulation (1/2)
    note:
      - If geom = 'circle', R0 is the radius.
      - If geom = 'square', R0 is the half side length.

**Subroutines** ``(vdet_calculate)``

.. code-block::
    :number-lines:

    subroutine vdet_calculate_current_1d
      Calculate and record the instantaneous momentum, probability current, probability density, and phase for a VD.
      Linear interpolation is used to obtain wavefunction quantities at the VD location.
    input
      this: virtual detector class object
      wavefn: Schrodinger 1D wavefunction class object
      x: spatial mesh array
      k: current time-step index
    output
      this % Krt: instantaneous momentum
      this % Jrt: instantaneous probability current
      this % rho: instantaneous probability density
      this % phase: instantaneous wavefunction phase

.. code-block::
    :number-lines:

    subroutine vdet_calculate_current_2d
      Calculate and record the instantaneous momentum, probability current, and phase for a VD.
      Bicubic interpolation is used to obtain wavefunction quantities at the VD location.
    input
      this: virtual detector class object
      wavefn: Schrodinger 2D wavefunction class object
      IDM: interpolation data matrix
      x, y: spatial mesh arrays
      nr: no. spatial grid points
      k: current time-step index
    output
      this % Krt, Jrt: instantaneous momentum and probability current
      this % phase: instantaneous wavefunction phase at the VD position

.. code-block::
    :number-lines:

    subroutine vdet_calculate_current_from_phase
      Calculate and record the instantaneous momentum and probability current for a virtual detector
      using the phase of the wavefunction.
    input
      this: virtual detector class object
      psi: wavefunction
      dr: spatial step size
      nr: no. spatial grid points
      k: current time-step index
    output
      this % Krt, Jrt: instantaneous momentum and probability current
    note
      - This routine is NOT recommended over vdet_calculate_current
        due to the phase-unwrapping problem.
      - This routine is deprecated.

**Subroutines** ``(edet)``

.. code-block::
    :number-lines:

    subroutine edet_detect
      Record electron information.
    input
      this: end detector class object
      electron: group of electron class objects
      sdim: dimensionality of the simulation (1/2)
    output
      data: recorded electron information

**Subroutines** ``(electron_dt_propagate)``

.. code-block::
    :number-lines:

    subroutine electron_dt_propagate_analytic_nd
      Push an nD electron trajectory using a forward-Euler Hamiltonian integrator.
    input
      this: electron class object
      pdot: analytically pre-computed force vector
      dt: temporal step size
    output
      updated electron position and momenta

.. code-block::
    :number-lines:

    subroutine electron_dt_propagate_numeric
      Push an electron using a forward-Euler Hamiltonian integrator.
      Numerically differentiates a potential energy function.
    input
      this: electron class object
      V: potential energy function
      x, y: spatial mesh arrays
      dr, dt: spatial/temporal step size
      nr: no. spatial grid points
    output
      updated electron position and momenta
    note
      - An electron grid (class grid_electron) should be used to deposit
        the potential energy function and electron coordinates.
