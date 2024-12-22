MODULES
=======
Revised 12/12/2024

Mathematics ``(math.f08)``
--------------------------

**Functions** ``(linspace)``

.. code-block::
    :number-lines:

    function linspace_n
      Construct a linearly-spaced array of a specified length.
    input
      a, b: lower/upper interval bounds
      n: no. points
    output
      vec: array of values between [a,b]

.. code-block::
    :number-lines:

    function linspace_d
      Construct a linearly-spaced array of a specified step size.
    input
      a, b: lower/upper interval bounds
      dx: step size
    output
      vec: array of values between [a,b]

**Functions** ``(zeros/ones/Kdelta)``

.. code-block::
    :number-lines:

    function zeros
      Return an n-dimensional array of zeros.
    input
      n: no. points
    output
      vec: array of zeros with length n

.. code-block::
    :number-lines:

    function ones
      Return an n-dimensional array of ones.
    input
      n: no. points
    output
      vec: array of ones with length n

.. code-block::
    :number-lines:

    function Kdelta
      Kronecker delta function.
    input
      n, m: integers
    output
      r: one if n == m, and zero otherwise.

**Functions** ``(rotM90)``

.. code-block::
    :number-lines:

    function rotM90X
      Rotate a square matrix 90 degrees counter-clockwise.
    input
      M: input matrix
    output
      V: rotated matrix

**Functions** ``(fliplr)``

.. code-block::
    :number-lines:

    function fliplr_xa
      Reverse-order the entries of an array.
    input
      arr: input array
    output
      f_arr: flipped array

.. code-block::
    :number-lines:

    function fliplr_xm
      Reverse-order the elements of a matrix.
    input
      M: input matrix
      ax: axis to flip over
    output
      f_M: flipped matrix
    note:
      - For ax = 0, both rows and columns are flipped.
      - For ax = 1 (ax = 2), columns (rows) will be flipped upside-down (left-right).

**Functions** ``(Identity/diag/inverse/Kproduct)``

.. code-block::
    :number-lines:

    function Identity
      Returns the n-by-n identity matrix.
    input
      n: square matrix dimension
    output
      A: Identity matrix

.. code-block::
    :number-lines:

    function diag_x
      Creates a square diagonal matrix from an array.
    input
      arr: input array of pre-allocated length
    output
      M: diagonal matrix

.. code-block::
    :number-lines:

    function inverse
      Compute the inverse of a real matrix.
    input
      M: square matrix
    output
      C: inverse of M
    note
      - This algorithm is based on Doolittle LU decomposition for Ax=b.
      - If M is singular, C is a matrix of NaNs.

.. code-block::
    :number-lines:

    function Kproduct
      Returns the Kronecker product of two matrices.
    input
      A, B: input matrices of pre-allocated dimensions
    output
      AB: Kronecker product of A and B

**Subroutines** ``(tridiag)``

.. code-block::
    :number-lines:

    subroutine tridiag_matmul_dmat_sp_xyz
      Sparse multiplication of a tri-diagonal with a diagonal matrix.
    input
      A: tri-diagonal matrix
      D: diagonal matrix
    output
      B: solution matrix A.D (tri-diagonal)
    note
      - Matrix A must be input as a 3n-2 vector in row-major dense-to-sparse ordering.
      - Matrix D must be input as a size-n vector.
      - Matrix B is output as a 3n-2 vector, like A.
      - Invoke the correct version through the interface tridiag_matmul_dmat.

.. code-block::
    :number-lines:

    subroutine tridiag_matmul_cvec_sp_xyz
      Sparse multiplication of a tri-diagonal matrix by a column vector.
    input
      A: tri-diagonal matrix
      x: column vector
    output
      b: solution vector A.x
    note
      - Matrix A must be input as a 3n-2 vector in row-major dense-to-sparse ordering.
      - The output b is a size-n vector representing the product.
      - Invoke the correct version through the interface tridiag_matmul_cvec.

.. code-block::
    :number-lines:

    subroutine tridiag_fbwd_subs_sp_xyz
      Solve the tri-diagonal matrix equation Ax=b using sparse forward-backward substitution.
    input
      A, b: tri-diagonal (sparse) matrix, vector
    output
      x: solution vector
    note
      - Matrix A must be input as a 3n-2 vector in row-major dense-to-sparse ordering.
      - Vector b is modified upon execution.
      - Invoke the correct version through the interface tridiag_fbwd_subs.

**Subroutines / Functions** ``(svdcmp)``

.. code-block::
    :number-lines:

    subroutine svdcmp
      Singular value decomposition.
      Given a matrix A = a(m,n), this routine computes its singular value decomposition, A = U.W.Vt.
      The matrix U replaces A on output. The diagonal matrix of singular values W is output
      as a vector w(n), and the matrix V (not the transpose Vt) is output as v(n,n).
    input
      a: original input matrix
      m, n: dimensions of a
      w, v: output placeholders
    output
      a: transformed matrix U
      w, v: vector of singular values diag(W) and matrix V
    note
      - This subroutine over-writes its inputs (a,w,v), and it depends on the pythag function.

.. code-block::
    :number-lines:

    function pythag
      Computes sqrt(a**2 + b**2) without destructive underflow or overflow.
    input
      a, b: values
    output
      pythag: sqrt(a**2 + b**2)
    note
      - Used mainly by the svdcmp subroutine.

**Subroutines / Functions** ``(fft)``

.. code-block::
    :number-lines:

    function fft_freq
      Construct a frequency array associated with a spatial/temporal variable.
    input
      n: no. points
      delta: sampling rate
      shift: shift the zero-frequency component to the center of the spectrum (boolean)
    output
      freq: array of frequency values
    note
      - The unshifted array has the following order.
         Term #:           Frequency:
         1 through n/2     positive [from 0 to Nyquist]
         n/2+1 through n   negative [-Nyquist to 0)
      - To obtain the angular frequency, use: 2*pi*freq.

.. code-block::
    :number-lines:

    subroutine fft_shift_nd
      Shift the zero-frequency component of a Fourier array to the center of the spectrum.
    input
      func: complex array
    output
      func: zero-shifted array

.. code-block::
    :number-lines:

    subroutine four1
      One-dimensional Fast Fourier Transform (FFT) routine.
    input
      data: 1D array to Fourier transform
      nn: no. points
      isign: operation; forward (+1) or inverse (-1) transform
    output
      data: Fourier-transformed array
    note
      - data must be a real array of length 2*nn representing alternating real/imaginary
        parts of a complex array.
      - nn must be a power of 2.
      - if isign = -1, the result is multiplied by nn.
      - To interface, invoke the wrapping functions fft_1d and ifft_1d.
      - Refer to: W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery,
        Numerical Recipes in Fortran 90 (Cambridge University Press, Cambridge, 2001).

.. code-block::
    :number-lines:

    subroutine four2
      Two-dimensional Fast Fourier Transform (FFT) routine.
    input
      data: 2D matrix to Fourier transform
      nn: integer array of dimensions; dim(data)
      isign: operation; forward (+1) or inverse (-1) transform
    output
      data: Fourier-transformed matrix
    note
      - data must be a real array of length 2*nn(1)*nn(2) representing alternating
        real/imaginary parts of a complex matrix.
      - All elements of nn must be a power of 2.
      - if isign = -1, the result is multiplied by nn(1)*nn(2).
      - To interface, invoke the wrapping functions fft_2d and ifft_2d.
      - Refer to: W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery,
        Numerical Recipes in Fortran 90 (Cambridge University Press, Cambridge, 2001).

.. code-block::
    :number-lines:

    function fft_nd
      n-dimensional Fast Fourier Transform.
    input
      func: matrix to Fourier transform
    output
      ft_func: Fourier-transformed matrix
    note:
      - Every dimension of func must be a power of 2.
      - The result is unshifted in frequency space.

.. code-block::
    :number-lines:

    function ifft_nd
      n-dimensional Inverse Fast Fourier Transform.
    input
      func: matrix to inverse Fourier transform
    output
      ift_func: inverse Fourier-transformed matrix
    note:
      - Every dimension of func must be a power of 2.
      - The result is unshifted in frequency space.
      - The factor of size(func) arising in four1/2 is divided out.

.. code-block::
    :number-lines:

    function real_stagger_complex_nd
      Staggers a complex array into alternating real/imaginary parts.
    input
      g: complex array
      n: no. points
    output
      f: real array of staggered values
    note
      - Odd indices: Re(g), Even: Im(g) elements.

.. code-block::
    :number-lines:

    function complex_stagger_real_nd
      Staggers elements of a real array into a complex array.
      Inverse operation performed by real_stagger_complex_nd.
    input
      f: real array
      n: half no. points; length(f) = 2*n
    output
      g: complex array of staggered values
    note
      - For a natural number j, g(j) = f(2*j-1) + i*f(2*j).

.. code-block::
    :number-lines:

    subroutine zero_pad_signal_x1d
      Zero-pads a 1D array, extending its length to a specified size.
      The corresponding time array is also extended but its step size is preserved.
    input
      x: time data
      f: signal data
      n: new size after padding
      dx: time step size
    output:
      x, f: zero-padded arrays
    note
      - The step size dx is recalculated after extending x, though it should not change.
      - This subroutine is primarily used if the signal is to be Fourier transformed
        but its length is not a power of 2 (see function next_pow2).

.. code-block::
    :number-lines:

    logical function is_pow2
      Test if an integer is a power of 2.
    input
      n: integer

.. code-block::
    :number-lines:

    integer function next_pow2
      Returns the nearest power of 2.
    input
      n: target integer

**Functions** ``(deriv/grad)``

.. code-block::
    :number-lines:

    subroutine deriv_x1d
      Compute the numerical derivative of an array using a five-point stencil.
    input
      f: input function
      dx: grid step size
    output
      df: derivative array

.. code-block::
    :number-lines:

    subroutine deriv_x2d
      Compute the numerical derivative of a matrix along a specified axis using a five-point stencil.
    input
      f: input matrix
      dr: grid step size
      ax: differentiation axis (index)
    output
      df: derivative matrix
    note:
      - For ax = 1 (ax = 2), the difference along rows (columns) will be calculated.

.. code-block::
    :number-lines:

    function grad_x1d
      Compute the one-dimensional gradient of a function.
    input
      f: input function
      del: grid step size
    output
      gradf: gradient of the input function

.. code-block::
    :number-lines:

    function grad_x2d
      Compute the gradient of a function along a specified axis.
    input
      f: input function
      del: grid step size
      ax: differentiation axis (index)
    output
      gradf: gradient of the input function
    note
      - For ax = 1 (ax = 2), the difference along rows (columns) will be calculated.

**Functions** ``(diff/simint/trapz/phase_unwrap)``

.. code-block::
    :number-lines:

    function diff
      Compute the discrete difference between adjacent array elements.
    input
      x: input array
    output
      dx: difference array
    note
      - The length of dx is size(x)-1.

.. code-block::
    :number-lines:

    function simint
      Cumulative numerical integration using a modified Simpson's Rule.
    input
      y: real array to be integrated
      y0: initial value
      dx: spatial grid step size
    output
      inty: integral of y(x)
    note
      - Refer to: L. V. Blake, U.S. NRL Memorandum Report 2231 (1971), titled:
        "A Modified Simpson's Rule and Fortran Subroutine for Cumulative Integration
        of a Function Defined by Data Points"

.. code-block::
    :number-lines:

    function trapz_x1d
      Numerical integration of an array using the trapezoidal formula.
    input
      f: array to integrate
    output
      s: total integral of f

.. code-block::
    :number-lines:

    function trapz_x2d_part
      Numerical integration of a 2D matrix using the trapezoidal formula.
      Integrates through a particular axis to produce an array.
    input
      f: matrix to integrate
      ax: integration axis (index)
    output
      s: numerical integral of f along ax
    note
      - For ax = 1 (ax = 2), the sum along rows (columns) will be calculated.
      - The chosen axis is integrated *through*, so for example if
      dim(f) = (nx,ny) and ax = 1, then dim(s) = nx.

.. code-block::
    :number-lines:

    function trapz_x3d_part
      Numerical integration of a 3D matrix using the trapezoidal formula.
      Integrates out a particular axis to produce a 2D matrix.
    input
      f: matrix to integrate
      ax: integration axis (index)
    output
      s: numerical integral of f along ax
    note
      - The chosen axis is integrated *out*, so for example if
      dim(f) = (nx,ny,nz) and ax = 2, then dim(s) = (nx,nz).

.. code-block::
    :number-lines:

    function trapz_x2d_full
      Numerical integration of a 2D matrix using the trapezoidal formula.
      Integrates over all rows and columns.
    input
      f: matrix to integrate
    output
      s: total integral of f
    note
      - The output can be multiplied by the product of step sizes, dr(2).

.. code-block::
    :number-lines:

    function trapz_x3d_full
      Numerical integration of a 3D matrix using the trapezoidal formula.
      Integrates over all rows and columns.
    input
      f: matrix to integrate
    output
      s: total integral of f
    note
      - The output can be multiplied by the product of step sizes, dr(3).

.. code-block::
    :number-lines:

    subroutine phase_unwrap_nd
      Unwrap 2-pi phase jumps arising from the arctangent function (atan2).
    input
      f: input phase array/matrix
    output
      f: phase-unwrapped array/matrix

**Subroutines / Functions** ``(bcuint/LegendrePoly/winHann/init_RNG)``

.. code-block::
    :number-lines:

    subroutine bcuint_x
      Bicubic interpolation within a Cartesian mesh.
    input
      y, y1, y2, y12: function, gradients, and cross derivative at the four grid
                      points of a rectangular cell (numbered ccw from lower left)
      x1l, x2l: lower-bound points on the coarse grid closest to the interpolation
                point in the x1- and x2-direction
      x1, x2: interpolation point coordinates
      dx: coarse grid step sizes
    output
      ansy: interpolated function value
      ansy1, ansy2: interpolated gradient values
    note
      - This routine performs the same task as bcuint_r_old, though it is slightly more optimized.

.. code-block::
    :number-lines:

    subroutine bcuint_r_old
      Bicubic interpolation within a Cartesian mesh. Deprecated version.
    input
      y, y1, y2, y12: function, gradients, and cross derivative at the four grid
                      points of a rectangular cell (numbered ccw from lower left)
      x1l, x1u,...: lower/upper coordinates in the x1- and x2-direction
      x1, x2: interpolation point coordinates
    output
      ansy: interpolated function value
      ansy1, ansy2: interpolated gradient values
    note
      - This routine calls bcucof for the interpolation coefficients.
      - Refer to: W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery,
        Numerical Recipes in Fortran 90 (Cambridge University Press, Cambridge, 2001).

.. code-block::
    :number-lines:

    subroutine bcucof
      Coefficients for bicubic interpolation.
    input
      y, y1, y2, y12: see description for bcuint_r_old
      d1, d2: grid cell length in the x1- and x2-direction
    output
      c: table of coefficients used by the routine bcuint_r_old for bicubic interpolation

.. code-block::
    :number-lines:

    recursive function LegendrePoly
      Evaluate the n-th degree Legendre polynomial at the point x.
    input
      n, x: polynomial degree (>= 0) and evaluation point
    output
      r: P_n(x)
    note
      - The Legendre polynomials are defined by the recursion relation:
      P_0(x) = 1.0, P_1(x) = x, and (n+1)*P_n+1(x) = (2n+1)*x*P_n(x) - n*P_n-1(x).

.. code-block::
    :number-lines:

    function LegendrePolySeq
      Generate a sequence of Legendre polynomials evaluated at the point x.
    input
      n, x: polynomial degree (>= 2) and evaluation point
    output
      s: n-dimensional array of Legendre polynomials, [P_0(x),...,P_n-1(x)]
    note
      - The first element of the sequence is P_0(x) = 1.0, and the last element
      is the (n-1)th degree Legendre polynomial.

.. code-block::
    :number-lines:

    pure recursive function factorial
      Evaluate n-factorial.
    input
      n: integer
    output
      r: n!
    note
      - Accurate for n <= 33.

.. code-block::
    :number-lines:

    subroutine cache_factorial
      Cache the first 33 factorials.
    input / output
      r: length(33) integer(16) array

.. code-block::
    :number-lines:

    function winHann
      Evaluates the Hanning window function.
    input
      t, tau: current/total time
    note
      - Used in eigenstate distillation.

.. code-block::
    :number-lines:

    subroutine init_RNG
      Initialize the (pseudo) Random Number Generator by querying /dev/urandom for seeds.

Quantum mechanics ``(quantum.f08)``
-----------------------------------

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

Electromagnetic field ``(emfm.f08)``
------------------------------------

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

Virtual detector ``(vdm.f08)``
------------------------------

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
      init: initialize a virtual detector and its data arrays.
      trigger(1/2): calculate and record the instantaneous momentum, probability current, probability density, and phase.
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
      apush(1/2): analytic propagation routine
      npush: numeric propagation routine

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

Rochester ``(rochester.f08)``
-----------------------------

**Functions** ``(Vsc)``

.. code-block::
    :number-lines:

    pure function Vsc
      Evaluates the soft-core Coulomb potential.
    input
      Z(2): charge numbers
      x(:): evaluation point
      s: screening parameter

.. code-block::
    :number-lines:

    pure function DVsc
      Evaluates the partial derivative of the soft-core Coulomb potential.
    input
      Z(2): charge numbers
      x(:): evaluation point
      s: screening parameter
      j: differentiation component

Optimization ``(optimize.f08)``
-------------------------------

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
