Documentation
=============
Module: Optimization ``(optimize.f08)``
---------------------------------------
Revised 9/14/2025

Note: Compiling ``optimize.f08`` with 64-bit floats (double precision) is recommended.

**Function interfaces**

Functions passed to optimization procedures must conform to one of the following templates.

.. code-block::
    :number-lines:

    ! scalar->scalar function
    function ssfunc(x)
        use prec, only: num
        implicit none
        real(num), intent(in) :: x
        real(num) :: ssfunc
    end function ssfunc

.. code-block::
    :number-lines:

    ! vector->scalar function
    function vsfunc(x)
        use prec, only: num
        implicit none
        real(num), intent(in) :: x(:)
        real(num) :: vsfunc
    end function vsfunc

.. code-block::
    :number-lines:

    ! vector->vector function
    function vvfunc(x)
        use prec, only: num
        implicit none
        real(num), intent(in) :: x(:)
        real(num) :: vvfunc(size(x))
    end function vvfunc

**Optimize-ND classes**

.. code-block::
    :number-lines:

    type OptimizeND_NelderMead
      Nelder-Mead N-dimensional minimization class
    variables
      func: objective function (vector->scalar)
      y, p: evolving objective function values & simplex vertices
      ftol: target fractional tolerance
      ndim: no. independent variables
      iter: iteration counter
      itmax: maximum no. iterations
      warn: toggle warning if iter exceeds itmax
      aux(n): auxiliary matrices; store anything you'd like
    procedures (core)
      create => create_NM: set the target function and object parameters
      minimize => amoeba: execute the Nelder-Mead minimization routine
      destroy => destroy_NM: object destructor
      reset => reset_NM: zero optimization variables and iter
    procedures (support)
      amotry => amotry_NM: simplex extrapolation tester
    note
      - The variable dimensions are:
          y, p: (ndim+1), (ndim+1,ndim)
      - The ndim+1 rows of p are size(ndim) vectors identifying the simplex vertices.
      - The elements of y equal the target function (func) evaluated at the ndim+1 vertices (rows) of p.
      - Before calling minimize(), the user must set the initial values of (y,p).
      - Upon completion, (y,p) will be ndim+1 new points all within ftol of a minimum,
        and iter will equal the number of function evaluations taken.

.. code-block::
    :number-lines:

    type OptimizeND_ConjugateGradient
      Conjugate-gradient N-dimensional minimization class
      based on the Fletcher-Reeves-Polak-Ribiere (FRPR) algorithm
    variables
      func, gfunc: objective function (vector->scalar) & its gradient (vector->vector)
      p: evolving independent vector
      y, xi: objective function value and gradient vector at p
      alpha: guess of line-minimization bracketing extent
      ftol: target fractional tolerance
      ndim: no. independent variables
      iter: iteration counter
      itmax: maximum no. iterations
      warn: toggle warning if iter exceeds itmax
    procedures (core)
      create => create_CG: set the target function, its gradient, and object parameters
      minimize => frprmn: execute the FRPR conjugate-gradient minimization routine
      destroy => destroy_CG: object destructor
      reset => reset_CG: zero optimization variables and iter
    procedures (support)
      linmin => linmin_CG: line-minimization subroutine
      mnbrak => mnbrak_CG: specialized minimum bracketing
      dbrent => dbrent_CG: Brent's derivative-based method
    note
      - The variable dimensions are:
          xi & p: (ndim)
      - Before calling minimize(), the user must set the initial values of (y,p).
      - Upon completion, (y,p) will be a new function value/point within ftol of
        a minimum, and iter will equal the number of iterations taken.

.. code-block::
    :number-lines:

    type OptimizeND_ParticleSwarm
      Particle-swarm N-dimensional minimization class
    variables
      func: objective function (vector->scalar)
      y, p: final objective function value and array of independent variables
      x, v: instantaneous particle positions & momenta
      xbest: best particle positions found thus far
      func_xbest: corresponding best objective function values found thus far
      w: the inertia weight constant (optional, default: 0.8)
      c(2): the cognitive & social coefficients (optional, default: both 0.1)
      ftol: target tolerance
      ndim: no. independent variables
      npart: no. particles in the swarm
      iter: iteration counter
      itmax: maximum no. iterations
      parallel: toggle parallel (OMP) advancement of particles (optional, default: false)
      warn: toggle warning if iter exceeds itmax
    procedures (core)
      create => create_PS: set the target function and object parameters
      minimize => nemo: execute the particle-swarm minimization routine
      destroy => destroy_PS: object destructor
      reset => reset_PS: zero optimization variables and iter
    procedures (support)
      span => flockspan_PS: calculate the characteristic size of the swarm
    note
      The variable dimensions are:
        p: (ndim)
        func_xbest: (npart)
        x, v, xbest: (npart,ndim)
      - Before calling minimize(), the user must set the initial values of (x,v).
        The nemo subroutine performs first calls to the target function.
      - Upon completion, (y,p) will be a new function value/point within ftol of
        a minimum, and iter will equal the number of iterations taken.

.. code-block::
    :number-lines:

    type OptimizeND_GuidedMonteCarlo
      Guided Monte Carlo N-dimensional minimization class
    variables
      func: objective function (vector->scalar)
      p: evolving independent vector
      y: objective function value
      mag_step, _delta: guided step and perturbation size amplitude for ea. coordinate
      step: normalized guide lengths that decrease with increasing frustration
      iroc: maximum no. times to "rock" the variables
      ndim: no. independent variables
      iter: iteration counter
      itmax: maximum no. iterations
      warn: toggle warning if iter exceeds itmax
    procedures (core)
      create => create_GMC: set the target function and object parameters
      minimize => gmcmn: execute the guided Monte Carlo minimization routine
      destroy => destroy_GMC: object destructor
      reset => zero optimization variables and iter
    note
      - The variable dimensions are:
          p, mag_step, mag_delta: (ndim)
      - Before calling minimize(), the user must set the initial values of (y,p).
      - This algorithm is due to:
          R. Delgoda and J. D. Pulfer, "A Guided Monte Carlo Search Algorithm for Global
          Optimization of Multidimensional Functions", J. Chem. Inf. Comput. Sci. vol. 38, pp. 1087-1095 (1998).
        (The terminology used throughout is based on this paper.)

**Subroutines** ``(initialize)``

.. code-block::
    :number-lines:

    subroutine create_NM
      Set the target function & parameters for a Nelder-Mead minimization task
    input
      this: Nelder-Mead class object
      func: user-defined vector->scalar function (vsfunc template)
      ndim: no. independent variables
      ftol: desired fractional tolerance
      itmax: maximum no. iterations
      warn: toggle warning if iter exceeds itmax
    note
      - Default values:
          ftol=1.d-6 if 0. is passed.
          itmax=5000 if -1 is passed.

.. code-block::
    :number-lines:

    subroutine create_CG
      Set the target function & parameters for an FRPR conjugate-gradient minimization task
    input
      this: FRPR conjugate-gradient class object
      func: user-defined vector->scalar function (vsfunc template)
      gfunc: the vector->vector gradient of func (vvfunc template)
      alpha: guess of line-minimization bracketing extent
            (could be on the order of the abscissa length-scale)
      ndim: no. independent variables
      ftol: desired fractional tolerance
      itmax: maximum no. iterations
      warn: toggle warning if iter exceeds itmax
    note
      - Default values:
          ftol=1.d-6 if 0. is passed.
          itmax=200 if -1 is passed.

.. code-block::
    :number-lines:

    subroutine create_PS
      Set the target function & parameters for a particle-swarm minimization task
    input
      this: particle-swarm class object
      func: user-defined vector->scalar function (vsfunc template)
      ndim: no. independent variables
      npart: no. particles in the swarm
      ftol: desired tolerance
      itmax: maximum no. iterations
      warn: toggle warning if iter exceeds itmax
      parallel: toggle parallel (OMP) advancement of particles (optional)
      w: the inertia weight constant (optional)
      c(2): the cognitive & social coefficients (optional)
    note
      - Default values:
          w=0.8 if nothing is passed in
          c=(/0.1,0.1/) if nothing is passed in

.. code-block::
    :number-lines:

    subroutine create_GMC
      Set the target function & parameters for a Guided Monte Carlo minimization task
    input
      this: Guided Monte Carlo class object
      func: user-defined vector->scalar function (vsfunc template)
      mag_step, _delta: guided step and perturbation size amplitude for ea. coordinate
      ndim: no. independent variables
      itmax: maximum no. iterations
      iroc: maximum no. times to "rock" the variables
      warn: toggle warning if iter exceeds itmax
    note
      - Default values:
          itmax=10,000 if -1 is passed
          itmax is ignored (run to completion) if -2 is passed
          iroc=3 if -1 is passed
      - The values of mag_step(i) and mag_delta(i) should approximately equal the
        characteristic scale-length associated with the coordinate p(i). For example,
        if p(i) is an azimuthal angle, the objective function containing cos(p(i)) terms,
        then mag_step(i)=2.*pi and mag_delta(i)=2.*pi/10 would be reasonable. (The latter
        quantity is an order of magnitude smaller as "delta" is a small perturbation to p(i)
        in the core Guided Monte Carlo procedure.)

**Procedures** ``(Nelder-Mead)``

.. code-block::
    :number-lines:

    subroutine amoeba
      Core subroutine for Nelder-Mead N-dimensional minimization
    input
      this: Nelder-Mead class object
    output
      y, p: function values & simplex vertices of a minimum
      iter: no. function evaluations performed
    note
      - Depends on the amotry_NM function.

.. code-block::
    :number-lines:

    function amotry_NM
      Support procedure for Nelder-Mead N-dimensional minimization
    description
      Extrapolates by a factor "fac" through the face of the simplex across from the
      high point, tries it, and replaces the high point if the new point is better.
    input, output
      See subroutine amoeba

**Procedures** ``(FRPR conjugate-gradient)``

.. code-block::
    :number-lines:

    subroutine frprmn
      Core subroutine for FRPR conjugate-gradient minimization
    input
      this: FRPR conjugate-gradient class object
    output
      y, p: function value & point of a minimum
      iter: no. function evaluations performed
    note
      - We implement the Polak-Ribiere variant for computing a conjugacy
        coefficient; the original Fletcher-Reeves version is commented out.

.. code-block::
    :number-lines:

    subroutine linmin_CG
      Support procedure for FRPR conjugate-gradient minimization
    description
      Given a point p and direction xi, moves and resets p to where the function
      func(p) takes on a minimum along the direction xi from p, and replaces xi by
      the actual vector displacement that p was moved. Also returns the function value
      at the returned location p. This is all accomplished by calling the routines
      mnbrak_CG and dbrent_CG.
    input, output
      See subroutine frprmn
    note
      - The invoked mnbrak_CG and dbrent_CG procedures are wrappers for the standard one-dimensional versions,
        but with functions f1dim and df1dim (defined in Ch. 10 of Ref. 1) substituted for this FRPR method.

.. code-block::
    :number-lines:

    subroutine mnbrak_CG
      Support procedure for FRPR conjugate-gradient minimization
    description
      Equivalent to the standard one-dimensional mnbrak subroutine (see its description)
      but specialized to evaluate func(p + X*xi) for this FRPR class object, where X=scalar.
    input, output
      See subroutine linmin_CG

.. code-block::
    :number-lines:

    function dbrent_CG
      Support procedure for FRPR conjugate-gradient minimization
    description
      Equivalent to the standard one-dimensional dbrent function (see its description)
      but specialized to evaluate func/gfunc(p + X*xi) for this FRPR class object, where X=scalar.
    input, output
      See subroutine linmin_CG

**Procedures** ``(Particle-swarm)``

.. code-block::
    :number-lines:

    subroutine nemo
      Core subroutine for particle-swarm N-dimensional minimization
    input
      this: particle-swarm class object
    output
      y, p: function value & point of a minimum
      iter: no. function evaluations performed
    note
      - Depends on the flockspan_PS function.

.. code-block::
    :number-lines:

    function flockspan_PS
      Support procedure for particle-swarm N-dimensional minimization
    description
      Computes the swarm's characteristic size: the smallest L2 distance between any two particles.
    input, output
      See subroutine nemo
    note
      - A better metric may be the *largest* L2 distance between any two particles...

**Procedures** ``(Guided Monte Carlo)``

.. code-block::
    :number-lines:

    subroutine gmcmn
      Core subroutine for Guided Monte Carlo N-dimensional minimization
    input
      this: Guided Monte Carlo class object
    output
      y, p: function value & point of a minimum
      iter: no. function evaluations performed

**Optimize-1D core procedures**

.. code-block::
    :number-lines:

    subroutine golden
      Golden Section Minimization (1D)
    description
      Given a function (func) and bracketing triplet of abscissas (ax,bx,cx),
      this routine performs a golden section search for the minimum, isolating
      it to a fractional precision of about ftol. The abscissa of the minimum
      is returned as xmin, and the minimum function value is returned as fmin.
      Parameters: R & C are the golden ratios.
    input
      func: function (scalar->scalar) to optimize
      ax, bx, cx: bracketing triplet of abscissas
      ftol: fractional precision, no smaller than ~sqrt(epsilon(1.d0))
    output
      xmin, fmin: abscissa & ordinate of minimum
    note
      - Use subroutine mnbrak to ensure (ax,bx,cx) is a bracketing triplet.

.. code-block::
    :number-lines:

    subroutine brent
      Brent's Minimization Method (1D)
    description
      Given a function (func) and bracketing triplet of abscissas (ax,bx,cx),
      this routine isolates the minimum to a fractional precision of about ftol
      using Brent's method. The abscissa of the minimum is returned as xmin, and
      the minimum function value is returned as fmin. Parameters: ITMAX is the max
      allowed number of iterations; CGOLD is 1-minus the Golden Ratio; and ZEPS protects
      against trying to achieve a fractional accuracy for a minimum that is exactly zero.
    input
      func: function (scalar->scalar) to optimize
      ax, bx, cx: bracketing triplet of abscissas
      ftol: fractional precision, no smaller than ~sqrt(epsilon(1.d0))
    output
      xmin, fmin: abscissa & ordinate of minimum
    note
      - Use subroutine mnbrak to ensure (ax,bx,cx) is a bracketing triplet.

.. code-block::
    :number-lines:

    function dbrent
      Brent's Derivative-based Minimization Method (1D)
    description
      Given a function & its derivative (func & dfunc), and a bracketing
      triplet of abscissas (ax,bx,cx), this routine isolates the minimum
      to a fractional precision of about ftol using a modification of Brent's
      method that uses derivatives. The abscissa of the minimum is returned as
      xmin, and the minimum function value is returned as dbrent.
    input
      func, dfunc: function to optimize & its derivative (both scalar->scalar)
      ax, bx, cx: bracketing triplet of abscissas
      ftol: fractional precision, no smaller than ~sqrt(epsilon(1.d0))
    output
      xmin, dbrent: abscissa & ordinate of minimum
    note
      - Use subroutine mnbrak to ensure (ax,bx,cx) is a bracketing triplet.

**Procedures** ``(1D support)``

.. code-block::
    :number-lines:

    subroutine mnbrak
      Minimization Bracketer (1D)
    description
      Given a function (func) and distinct initial points (ax,bx), this routine searches
      in the downhill direction and returns new points (ax,bx,cx) that bracket a minimum
      of the function. Also returned are the function values at the 3 points, (fa,fb,fc).
      Parameters: GOLD is the default ratio by which successive intervals are magnified;
      GLIMIT is the maximum magnification allowed for a parabolic-fit step.
    input
      func: target function (scalar->scalar)
      ax, bx: search domain min & max
    output
      ax, bx, cx: bracketing triplet, ordered: ax<bx<cx or ax>bx>cx.
      fa, fb, fc: function evaluated at (ax,bx,cx), ordered: fb<fa and fb<fc.
    note
      - The input order of ax & bx does not matter.
      - Upon output, ax<bx<cx or ax>bx>cx, and fb<fa and fb<fc always.

**Procedures** ``(Search)``

.. code-block::
    :number-lines:

    recursive function BinarySearch
      Binary search function
    description
      Finds the index of element x in array arr,
      or the index of the nearest value less than x.
    input
      arr: array to search
      ia, ib: indices to search between [start, end]
      x: target element
    output
      ir: index result, arr(ir) == x or arr(ir) < x.
    note
      - The array (arr) must be ordered (monotonically increasing)
      - Beware of edge cases.
          If x <= arr(1), then ir = 1.
          If x >= arr(size(arr)), then ir = size(arr).
        It is advised to check these conditions externally, before invoking BinarySearch.

**Procedures** ``(misc)``

.. code-block::
    :number-lines:

    function dfridr_ss (or _vs)
      Compute the numerical derivative of a scalar->scalar (or vector->scalar)
      function at a point using Ridders' method of polynomial extrapolation.
    input
      func: function to differentiate
      n: component to differentiate (_vs only)
      x: evaluation point
      h: estimated initial step-size
      err: returned error estimate
    note
      - The step-size (h) need not be small; it should be an increment in x over
        which func changes substantially. This procedure iteratively reduces h.
      - Parameters: CON is the step-size reduction factor per iteration; NTAB is
        the maximum tableau size; and the procedure returns when the error is
        SAFE worse than the best so far.
      - You can invoke me through the module interface dfridr.
