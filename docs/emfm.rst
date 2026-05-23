Documentation
=============
Module: Electromagnetic Field ``(emfm.f08)``
--------------------------------------------
Revised 9/14/2025

**Classes**

.. code-block::
    :number-lines:

    type emf
      Electromagnetic Field class
    parameters
      profile: temporal envelope (string)
      E0: field amplitude
      omg0: central frequency
      eps: ellipticity value
      ch1: linear chirp coefficient
      Ncyc_rf, _pl: no. cycles rise/fall & plateau
      CEP: carrier-envelope phase
      t_on: start time
      t_off: end time
      T0: period
      Tp: total duration
      Tpk: gaussian peak time
      Tfwhm: gaussian intensity FWHM duration
      it_on: start time index
      it_off: end time index
    variables
      Ex, Ey: electric field components
      Ax, Ay: vector potential components
      Cx, Cy: excursion components
    procedures
      init_trapz: initialize field with trapezoidal pulse envelope
      init_sine2: initialize field with sine-squared pulse envelope
      init_gauss: initialize field with gaussian pulse envelope
      init_gauss_l: initialize field with linearly-ramped gaussian pulse envelope
    note
      - The available profiles are:
          trapz-N: trapezoidal, N-cycle rise/fall; N of type float
          sine2: sine-squared
          gauss: gaussian
          gauss-l: gaussian with linearly-ramped wings
      - When creating a trapz-N pulse, the user must define:
          E0, omg0, eps, ch1, CEP, t_on, Ncyc_rf, Ncyc_pl
      - When creating a sine2 pulse, the user must define:
          E0, omg0, eps, ch1, CEP, t_on, Tp
      - When creating a gauss pulse, the user must define:
          E0, omg0, eps, ch1, CEP, Tpk, Tfwhm

**Subroutines** ``(emf)``

.. code-block::
    :number-lines:

    subroutine emf_trapezoidal_pulse
      Initialize EM-field with a 2-cycle turn-on/off trapezoidal temporal profile.
    input
      this: emf class object
      t: time array
    output
      this % Ex, Ey: field components
      this % Ax, Ay: vector potential components
      this % Cx, Cy: excursion components
    note
      - The amplitude is normalized by the ellipticity value.
      - If ch1 is non-zero, the pulse will be linearly chirped.

.. code-block::
    :number-lines:

    subroutine emf_sine_squared_pulse
      Initialize EM-field with a sine-squared temporal profile.
    input, output
      see emf_trapezoidal_pulse
    note
      see emf_trapezoidal_pulse

.. code-block::
    :number-lines:

    subroutine emf_gaussian_pulse / emf_gaussian_l_pulse
      Initialize EM-field with a gaussian temporal profile.
    input, output
      see emf_trapezoidal_pulse
    note
      - The gaussian_l envelope is linearly ramped between [2,3]*w0t where w0t = 1/e^2 radius.
        The field parameters must accommodate these ramps in the simulation time domain.
