Documentation
=============
Module: QED Monte-Carlo Particle ``(qed_mcpm.F08)``
---------------------------------------------------
Revised 2/28/2026

.. raw:: html

  <style> .red {color:#ff0000; font-weight:normal;} </style>
.. role:: red

**References**

[:red:`1`] C. P. Ridgers, J. G. Kirk, R. Duclous, T. G. Blackburn, C. S. Brady, K. Bennett, T. D. Arber, and A. R. Bell, \"Modelling gamma-ray photon emission and pair production in high-intensity laser–matter interactions,\" J. Comp. Phys. **260**, 273-285 (2014).
`DOI 10.1016/j.jcp.2013.12.007 <https://doi.org/10.1016/j.jcp.2013.12.007>`_

[:red:`2`] J. G. Kirk, A. R. Bell, and I. Arka, \"Pair production in counter-propagating laser beams,\" Plasma Phys. Control. Fusion **51**, 085008 (2009). `DOI 10.1088/0741-3335/51/8/085008 <https://doi.org/10.1088/0741-3335/51/8/085008>`_

**Compilation Notes**

.. line-block::
  The following compilation flags are available.
    ``-DNOWARN``: Suppress warning messages during runtime.
    ``-DUNSAFE``: Do not perform some logical error checks.
  These flags are omitted by default, which is the recommended operation unless you are certain of, and have rigorously tested, your implementation of the module.
  The use of the ``-DUNSAFE`` flag may result in segmentation faults and/or erroneous particle behavior.

**Variables**  ``(Units)`` | ``subroutine set_qedmcpm_units``.

.. code-block::
    :number-lines:

    lam0_cm: Normalization length [cm] (defined).
    k0_icm: Wavenumber [1/cm] (derived).
    w0_Hz: Frequency [Hz] (derived).
    hbar_: Reduced Planck constant [sim.u.] (derived).
    q_: Electric charge [sim.u.] (derived).
    F_: Characteristic electric/magnetic field strength [sim.u.] (derived).

To convert from simulation units [sim.u.] to the centimeter-gram-second [CGS] system, multiply by the corresponding factor below.

.. code-block::
    :number-lines:

    [length] = 1/k0_icm
    [mass] = m_e
    [time] = 1/w0_Hz
    [charge] = |e|
    [speed] = w0_Hz/k0_icm = c
    [momentum] = m_e*c
    [energy] = m_e*(c^2)
    [field] = w0_Hz*m_e*c/|e|

Notes:

- In expressions, ``m_e``, ``c``, and ``|e|`` are unity, while ℏ becomes the ``hbar_`` variable. For example, the QED critical field ``Ec [CGS] = (m_e^2)*(c^3)/(|e|*ℏ)`` is simply ``1/hbar_`` in ``[sim.u.]``.

- ``hbar_`` and ``q_`` are related by the fine-structure constant.

**Variables**  ``(Parameters)`` | ``subroutine set_qedmcpm_params``.

.. code-block::
    :number-lines:

    Ep_min: Threshold photon energy for emission [sim.u.] (default: 1.0).
    dt_max: Maximum time-step [sim.u.] for both adaptive and fixed stepping (default: 1.0).
    dt_scl: Adaptive time-step scale factor, in (0,1] (default: 1.d-2).
    photon_dynamics: Toggle photon propagation and e-/e+ pair production (default: .true.).

Notes:

- ``Ep_min = 1.0`` corresponds to a minimum photon energy of ``m_e*(c^2)`` ``(~511 keV)``.
- ``dt_max = 1.0`` corresponds to the (inverse) maximum synchrotron emission rate in a field of strength ``pi/(5.24*sqrt(3)*alpha_f) [sim.u.]``, or ``E ~ 1.2e-4*Ec``.

**Classes**

.. code-block::
    :number-lines:

    type particle
      Particle Data class
    variables
      m, q: particle mass, electric charge
      x(0:3), p(0:3): 4-position, 4-momentum vector
      dt: particle-push time-step
      qp: quantum parameter
      x0_, p0_: time, energy at creation
      stat: status integer
      id: identifier number (auxiliary)
    procedures
      push => push_particle_rk4, push_photon: propagate the particle forth by one time-step dt
    note
      - The particle status codes are as follows.
          stat = 1: Propagate and process normally.
          stat = 0: Halt propagation and physical processes.
          stat = -1: Annihilate particle. (Deprecated)
      - The variables x0_ and p0_ are set in the new_particle constructor function.

.. code-block::
    :number-lines:

    type pnode_t
      Particle Node class (linked-list member)
    variables
      par: object of type(particle)
      prev, next: pointers of type(pnode_t)

.. code-block::
    :number-lines:

    type plist_t
      Particle List class (linked list)
    variables
      head, tail: first & last elements of the linked list, pointers of type(pnode_t)
      N: total number of particles in the list
    procedures (core)
      create => create_particle: append a Particle Node to the list
      delete => delete_particle: delete a Particle Node from the list
      reset => reset_list: reset the list by deallocating all Particle Nodes
      push => list_push: propagate all particles in the list forth by one time-step
      emit => list_emit: process emission events in one time-step for all particles in the list
    procedures (support)
      count => count_particles: count the number of Particle Node objects in the list
      has_stat => exists_particle_wstat: check if the list contains a particle of a given status code
      get_id => get_particle_ids: returns an integer array of Particle id numbers
      get_charge => get_particle_charges: returns an array of Particle charges
      get_energy => get_particle_energies: returns an array of Particle energies
      get_qparam => get_particle_quantum_parameters: returns an array of Particle quantum parameters
      get_ivalue => get_particle_initial_quantities: returns an array of Particle creation times/energies
      get_fvalue => get_particle_final_quantities: returns an array of Particle final positions/momenta
    note
      - Do not mix massive (m>0) with massless (m=0) particles in a single plist_t instance.
        Putting massive particles of different mass (m) and charge (q) values is OK.
        E.g., putting (m,q)=(1,-1) (an electron) and (m,q)=(207,+1) (a muon+) in a single plist_t
        object is OK, but (m,q)=(0,0) (a photon) with (m,q)=(1,-1) will break the core procedures.

.. code-block::
    :number-lines:

    type diag_t
      Diagnostic data structure
    variables
      (min/max)_qp(1/2): the rank-global minimum/maximum lepton/photon (index 1/2)
        quantum parameter encountered during evaluation

      (min/max)_prob(1/2): the rank-global minimum/maximum photon/lepton-pair (index 1/2)
        emission probabilities encountered during evaluation

.. code-block::
    :number-lines:

    type field
      Electric & Magnetic field components class
    variables
      E(x,y,z)/B(x,y,z): electric/magnetic field component functions
    procedures (core)
      set => set_field_components: set the field pointer variables
      unset => unset_field_components: nullify the field pointer variables
    procedures (support)
      E => get_field_electric: get vector electric field at (r,t)
      B => get_field_magnetic: get vector magnetic field at (r,t)
    note
      - E(x,y,z)/B(x,y,z) must conform to the emfunc interface

.. code-block::
    :number-lines:

    type force
      Force field function class
    variables
      eval: vector force-field function
    procedures
      set => set_force_function: set the force-field pointer variable
      unset => unset_force_function: nullify the force-field pointer variable
    note
      - eval must conform to the frfunc interface

.. code-block::
    :number-lines:

    type qmc_table_t
      Quantum Monte-Carlo lookup tables
    variables
      path: path to lookup data
    variables (lookup arrays)
      eta(0,2): quantum parameter, electron (e-) or positron (e+)
      chi(1,3): quantum parameter, photon (y)
      frac1: fractional energy share, e-/e+ production
      h2: synchrotron function, y emission rate
      g3: emissivity function, e-/e+ production rate
    variables (lookup matrices)
      chi0: element (i,j) is chi(i,j) corresponding to eta0(i).
      Py0: element (i,j) is the cumulative probability of emitting a chi0(i,j) photon by an eta0(i) fermion.
      Pf1: element (i,j) is the cumulative probability of a chi1(i) photon respectively giving
        fractions frac1(j) and 1.0-frac1(j) of its energy to the e- and e+ during production.
    procedures
      load => load_qmc_tables: load lookup data from binary files into class variables
    note
      - Class variables ending in the same number (e.g., eta0, chi0, and Py0) belong together.
      - Refer to: [1] C. P. Ridgers et al., "Modelling gamma-ray photon emission and pair production
        in high-intensity laser–matter interactions", J. Comp. Phys. 260, 273 (2014).
      - The default lookup tables were obtained from a script that reproduces those in Ref. 1, Appendix D.

**Function interfaces**

Functions defining the electric-magnetic and force fields experienced by particles must respectively conform to the ``emfunc`` and ``frfunc`` templates outlined below.

Use ``emfunc`` to define a single field component (``E{x,y,z}``, ``B{x,y,z}``) as a function of space ``r(3)`` and time ``t``.

.. code-block::
    :number-lines:

    ! field function template
    function emfunc(r,t)
        use prec, only: num
        implicit none
        real(num), intent(in) :: r(3), t
        real(num) :: emfunc
    end function emfunc

Use ``frfunc`` to define a vector force field acting on a particle (mass ``m`` / charge ``q``) as a function of space ``r(3)``, momentum ``p(3)``, time ``t``, and electric-magnetic components contained within a ``field`` class object.

.. code-block::
    :number-lines:

    ! force function template
    function frfunc(m,q,r,p,t,F)
        use prec, only: num
        import field
        implicit none
        real, intent(in) :: m, q
        real(num), intent(in) :: r(3), p(3), t
        class(field), intent(in) :: F
        real(num) :: frfunc(3)
    end function frfunc

**Procedures** ``(set_qedmcpm)``

.. code-block::
    :number-lines:

    subroutine set_qedmcpm_units
      Set unit quantities (length, time, field, etc.) for a qed_mcpm calculation.
    input
      lam0_cm_: defining scale-length, in centimeters

.. code-block::
    :number-lines:

    subroutine set_qedmcpm_params
      Set parameters for a qed_mcpm calculation.
    input
      Ep_min_: threshold photon energy for emission (optional)
      dt_max_: maximum time-step (optional)
      dt_scl_: adaptive time-step scale factor in (0,1] (optional)
      photon_dynamics_: toggle photon propagation and e-/e+ pair production (optional)
    note
      - The module default values are listed above.

**Procedures** ``(push)``

.. code-block::
    :number-lines:

    subroutine push_particle_rk4
      Propagate a massive particle forth by one time-step using the 4th-order Runge-Kutta method.
    input
      par: Particle Data class object
      G: force class object
      F: field class object
    note
      - The time-step used is par%dt, which must be set externally.

.. code-block::
    :number-lines:

    subroutine push_particle_vay
      Propagate a massive particle forth by one time-step using J.-L. Vay's relativistic method.
    input
      par: Particle Data class object
      G: force class object
      F: field class object
    note
      - The time-step used is par%dt, which must be set externally.
      - This is a leapfrog scheme in which momentum (p) is half-advanced by par%dt/2 over position (x).
      - See J.-L. Vay, "Simulation of beams or plasmas crossing at relativistic velocity,"
        Phys. Plasmas 15, 056701 (2008). DOI https://doi.org/10.1063/1.2837054

.. code-block::
    :number-lines:

    subroutine push_photon
      Propagate a massless particle (photon) forth by one time-step.
    input
      par: Particle Data class object
    note
      - The time-step used is par%dt, which must be set externally.

**Procedures** ``(plist_t-create/delete)``

.. code-block::
    :number-lines:

    subroutine create_particle
      Instantiate and append a Particle Node to the list.
    input
      lst: Particle List class object
      par: Particle Data class object
    note
      - Usage example:
          type(plist_t) :: MyList
          MyList%create(new_particle(m=1., q=-1., id=7, x=[...], p=[...]))
        creates a particle of mass 1, charge -1, and ID tag 7 at 4-position x with 3-momentum p.
      - Refer to the new_particle constructor function.

.. code-block::
    :number-lines:

    subroutine delete_particle
      Delete a Particle Node from the list.
    input
      lst: Particle List class object
      node: Particle Node class object
    note
      - If the list is empty, the procedure simply returns (without notice).
      - If `node` is in the list, then it will be nullified (disassociated) upon return.
        Thus, if using delete_particle in a pointer-loop, it is critical to store the address of `node%next`.

.. code-block::

    subroutine reset_list
      Reset the list by deallocating all Particle Nodes.
    input
      lst: Particle List class object
    note
      - The particle count is also reset to zero.

.. code-block::
    :number-lines:

    function new_particle
      Particle Data constructor function.
    input
      m: particle mass
      q: electric charge
      id: identifier number (auxiliary)
      x(0:3): 4-position vector
      p(3), optional: 3-momentum vector
      E, optional: energy
      N, optional: direction 3-vector
      dt, optional: particle-push time-step
      stat, optional: status integer
    output
      par: a new Particle Data class object
    note
      - The 3-momentum (p) *or* relativistic energy (E) and direction-vector (N) must be set.
      - The direction 3-vector (N) need not be pre-normalized, but it cannot have zero norm.
      - The relativistic energy (E) must be greater than or equal to the rest-mass (m).
      - If E=m, then p(1:3) = 0.0, i.e., the particle will be initialized at rest.
      - If the particle-push time-step (dt) is omitted, it defaults to 1.d-3.
      - If the status integer (stat) is omitted, it defaults to 1 (normal).

**Procedures** ``(plist_t-push/emit)``

.. code-block::
    :number-lines:

    subroutine list_push
      Propagate all particles (nodes) in a Particle List forth by one time-step.
    input
      lst: Particle List class object
      t_max, optional: maximum particle time
      enrg_min, optional: minimum particle energy
      dt, optional: particle-push time-step
      G, optional: force class object
      F, optional: field class object
    note
      - If the particle-push time-step (dt) is omitted, then a field class object (F)
        must be provided in order for the adaptive time-step calculator to function.
      - For massive particles (m /= 0.0), force (G) and field (F) objects are required.
      - In this subroutine, the particle-push time-step is set (whether passed manually,
        or calculated through a dt_QED_ function) and reused in the subroutine list_emit.
      - A particle will cease to propagate / participate in processes (stat = 0) if any
        of the following conditions are met:
          (1) t_max is provided, and a particle's local time is greater than or equal to it.
          (2) enrg_min is provided, and a particle's energy is less than or equal to it.
          (3) dt is omitted, and the adaptive time-step calculator returns Infinity (the zero-field condition).

.. code-block::
    :number-lines:

    subroutine list_emit
      Process emission events for all particles (nodes) in a Particle List.
    input
      emit_lst: Particle List class object of the emitting species
      recv_lst: Particle List class object of the generated species
      F: field class object
      L: qmc_table_t class object
    note
      - The rate calculation is based on the per-particle time-step set in the subroutine list_push.
      - If emit_lst = photons, then recv_lst = leptons (electrons/positrons), and vice-versa.
      - Pair-production is energy-conserving, *not* momentum-conserving.
        (One cannot satisfy both simultaneously without a third body.)
      - If a photon energy is below 2.0, then pair production is energetically impossible.
        The photon will cease to propagate / participate in processes (stat = 0).
      - This function uses the relativistic small-emission-cone-angle approximation,
        i.e., the photon is emitted parallel to the lepton direction of propagation.

**Procedures** ``(plist_t-diagnostics)``

.. code-block::
    :number-lines:

    integer function count_particles
      Returns the number of Particle Nodes in the linked list.
    input
      lst: Particle List class object
    note
      - This function is not really needed; the particle count in a linked list
        is tracked and updated automatically with the (integer) class variable N.

.. code-block::
    :number-lines:

    function exists_particle_wstat
      Check if the linked list contains a particle of a given status code.
    input
      lst: Particle List class object
      stat: status integer
    output
      ans: true (false) if the list does (does not) contain a particle with the status `stat`

.. code-block::
    :number-lines:

    function get_particle_ids
      Returns an integer array of Particle id numbers in the linked list.
    input
      lst: Particle List class object
    output
      ans: integer array of Particle id numbers
    note
      - If the list is empty, `ans` will be a size-0 array.

.. code-block::
    :number-lines:

    function get_particle_charges
      Returns an array of Particle charges in the linked list.
    input
      lst: Particle List class object
    output
      ans: array of Particle charges
    note
      - If the list is empty, `ans` will be a size-0 array.

.. code-block::
    :number-lines:

    function get_particle_energies
      Returns an array of Particle energies in the linked list.
    input
      lst: Particle List class object
    output
      ans: array of Particle energies
    note
      - If the list is empty, `ans` will be a size-0 array.

.. code-block::
    :number-lines:

    function get_particle_quantum_parameters
      Returns an array of Particle quantum parameters in the linked list.
    input
      lst: Particle List class object
    output
      ans: array of Particle quantum parameters
    note
      - If the list is empty, `ans` will be a size-0 array.

.. code-block::
    :number-lines:

    function get_particle_initial_quantities
      Returns an array of Particle creation times/energies in the linked list.
    input
      lst: Particle List class object
      which: 'time'/'energy' for the creation time/energy
    output
      ans: array of Particle creation times/energies
    note
      - If the list is empty, `ans` will be a size-0 array.
      - Accepted aliases for time and energy are 't'/'T' and 'e'/'E'.

.. code-block::
    :number-lines:

    function get_particle_final_quantities
      Returns an array of Particle final positions/momenta in the linked list.
    input
      lst: Particle List class object
      which: 'position'/'momentum' for the final position/momentum
      ic: component index, in {0,1,2,3}
    output
      ans: array of Particle final positions/momenta
    note
      - If the list is empty, `ans` will be a size-0 array.
      - Accepted aliases for position and momentum are 'x'/'X' and 'p'/'P'.

**Procedures** ``(field)``

.. code-block::
    :number-lines:

    subroutine set_field_components
      Set pointers to user-defined field component functions.
    input
      F: field class object
      E(x,y,z)/B(x,y,z): electric/magnetic field component functions, conforming to the emfunc interface

.. code-block::
    :number-lines:

    subroutine unset_field_components
      Nullify pointers to field component functions.
    input
      F: field class object

.. code-block::
    :number-lines:

    function get_field_electric
      Get vector electric field at (r,t).
    input
      F: field class object
      r(3): position vector
      t: current time
    output
      E(3): electric field array

.. code-block::
    :number-lines:

    function get_field_magnetic
      Get vector magnetic field at (r,t).
    input
      F: field class object
      r(3): position vector
      t: current time
    output
      B(3): magnetic field array

**Procedures** ``(force)``

.. code-block::
    :number-lines:

    subroutine set_force_function
      Set pointer to user-defined force-field function.
    input
      G: force class object
      eval, optional: force-field function, conforming to the frfunc interface
    note
      - If eval is omitted, it defaults to the Force_Lorentz function.

.. code-block::
    :number-lines:

    subroutine unset_force_function
      Nullify pointer to force-field function.
    input
      G: force class object

**Procedures** ``(qmc_table_t)``

.. code-block::
    :number-lines:

    subroutine load_qmc_tables
      Load lookup data from binary files.
    input
      L: qmc_table_t class object
      path: path to lookup data

**Procedures** ``(quantum parameter)``

.. code-block::
    :number-lines:

    function Calc_QP_eta
      Calculate the electron/positron (e-/e+) quantum parameter.
    input
      r(3), p(3): position, momentum vector
      t: current time
      F: field class object
    note
      - The function implicitly uses a mass of unity, i.e., it
        is the quantum parameter of electrons/positrons only.

.. code-block::
    :number-lines:

    function Calc_QP_chi
      Calculate the photon quantum parameter.
    input
      r(3), p(3): position, momentum vector
      t: current time
      F: field class object
    note
      - The function includes a factor of 1/2, per the convention of Ref. 1.

**Procedures** ``(QED rates)``

.. code-block::
    :number-lines:

    function Calc_Rate_Synchrotron
      Calculate the rate of synchrotron emission by an electron/positron (e-/e+)
      via interpolation of qmc_table_t-type lookup data.
    input
      eta, enrg: e-/e+ quantum parameter, energy
      L: qmc_table_t class object
    output
      rate: synchrotron rate, with the units of inverse-time [sim.u.].
    note
      - The function implicitly uses a mass of unity,
        i.e., it is for electrons/positrons only.

.. code-block::
    :number-lines:

    function Calc_Rate_PairProduction
      Calculate the rate of pair-production by a photon
      via interpolation of qmc_table_t-type lookup data.
    input
      chi, enrg: photon quantum parameter, energy
      L: qmc_table_t class object
    output
      rate: pair-production rate, with the units of inverse-time [sim.u.].

**Procedures** ``(new-photon/pair)``

.. code-block::
    :number-lines:

    function Calc_newPhoton_Energy
      Calculate the energy of a photon newly-emitted by an electron/positron (e-/e+)
      via interpolation of qmc_table_t-type lookup data.
    input
      eta: quantum parameter (of the emitting e-/e+)
      par: Particle Data class object (of the emitting e-/e+)
      F: field class object
      L: qmc_table_t class object
    output
      enrg: new-photon energy
    note
      - This function uses the relativistic small-emission-cone-angle approximation,
        i.e., the photon is emitted parallel to the lepton direction of propagation.
      - Be mindful of edge-cases in lookup data interpolation,
        for quantum parameters that are above/below range.

.. code-block::
    :number-lines:

    function Calc_newPair_Fraction
      Calculate the fractional energy split of a photon into an electron/positron (e-/e+) pair
      via interpolation of qmc_table_t-type lookup data.
    input
      chi: quantum parameter (of the annihilating photon)
      L: qmc_table_t class object
    output
      frac: fractional energy split in [0,1]
    note
      - The e-/e+ receives (randomly) either frac or (1-frac) of the parent photon energy.
      - Low-energy photons are likely to split their energy evenly (frac = 1/2).
      - Be mindful of edge-cases in lookup data interpolation,
        for quantum parameters that are above/below range.

**Procedures** ``(utilities)``

.. code-block::
    :number-lines:

    function Norm3
      Return the magnitude of a 3-vector.
    input
      a(3): input 3-vector
    output
      Norm3: magnitude, sqrt(a.dot.a)

.. code-block::
    :number-lines:

    function Norm4
      Return the squared magnitude of a 4-vector.
    input
      a(0:3): input 4-vector
      sign: sign of the time component (integer, only +/-1)
    output
      Norm4: squared magnitude, sign*(time^2-space^2)
    note
      - Be mindful that this is the squared norm.

.. code-block::
    :number-lines:

    function Lorentz
      Evaluate the Lorentz (relativistic) factor for a particle.
    input
      p(3): momentum vector
      m: particle mass
    output
      Lorentz: relativistic factor

.. code-block::
    :number-lines:

    function Force_Lorentz
      Evaluate the relativistic Lorentz force for a particle at (r,t).
    input
      m, q: particle mass, electric charge
      r(3), p(3): position, momentum vector
      t: current time
      F: field class object
    output
      force(3): force array

.. code-block::
    :number-lines:

    function maxrate_QED_SYNCH
      Calculate the maximum synchrotron rate expected of an
      electron/positron (e-/e+) in a given electric-magnetic field.
    input
      E(3), B(3): electric, magnetic field array
    output
      maxrate_QED_SYNCH: maximum synchrotron rate, with the units of frequency [sim.u.]
    note
      - See the dt_QED_SYNCH function.
      - If Max{|E|,|B|} = 0.0, then maxrate_QED_SYNCH = 0.

.. code-block::
    :number-lines:

    function maxrate_QED_PPROD
      Calculate the maximum pair-production rate expected
      of a photon in a given electric-magnetic field.
    input
      E(3), B(3): electric, magnetic field array
    output
      maxrate_QED_PPROD: maximum pair-production rate, with the units of frequency [sim.u.]
    note
      - See the dt_QED_PPROD function.
      - If Max{|E|,|B|} = 0.0, then maxrate_QED_PPROD = 0.

.. code-block::
    :number-lines:

    function dt_QED_SYNCH
      Calculate the inverse of the maximum synchrotron rate expected
      of an electron/positron (e-/e+) in a given electric-magnetic field.
    input
      E(3), B(3): electric, magnetic field array
    output
      dt_QED_SYNCH: inverse synchrotron rate, with the units of time [sim.u.]
    note
      - Inverse of the maxrate_QED_SYNCH function.
      - Photon emission calculations should employ time-steps less than dt_QED_SYNCH,
        in order for the point-like formation length approximation to hold well.
      - If Max{|E|,|B|} = 0.0, then dt_QED_SYNCH = Infinity.
        (A fail-safe for this condition should be implemented when using dt_QED_SYNCH.)

.. code-block::
    :number-lines:

    function dt_QED_PPROD
      Calculate the inverse of the maximum pair-production rate
      expected of a photon in a given electric-magnetic field.
    input
      E(3), B(3): electric, magnetic field array
    output
      dt_QED_PPROD: inverse pair-production rate, with the units of time [sim.u.]
    note
      - Inverse of the maxrate_QED_PPROD function.
      - Pair-production calculations should employ time-steps less than dt_QED_PPROD,
        in order for the point-like formation length approximation to hold well.
      - If Max{|E|,|B|} = 0.0, then dt_QED_SYNCH = Infinity.
        (A fail-safe for this condition should be implemented when using dt_QED_PPROD.)
      - For all (E,B), dt_QED_PPROD > dt_QED_SYNCH.
        It is therefore practical to bound both processes temporally by dt_QED_SYNCH.

**Procedures** ``(error/warning)``

.. code-block::
    :number-lines:

    subroutine ERROR
      Display an error message and terminate the program with exit code 1.
    input
      message: the error message
      in, optional: the name of the calling subroutine, function, or file
      fix, optional: the recommended action to fix the error

.. code-block::
    :number-lines:

    subroutine WARNING
      Display a warning message and resume program execution.
    input
      message: the warning message
      in, optional: the name of the calling subroutine, function, or file
      effect, optional: the consequence of the warning
