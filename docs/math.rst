Documentation
=============
Module: Mathematics ``(math.f08)``
----------------------------------
Revised 9/14/2025

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

**Operators**

.. code-block::
    :number-lines:

    operator (.cross.)
      Cross-product between two vectors a(3) and b(3), c(3) = (a.cross.b).
    note
      - Use parentheses to ensure the correct order of operations.

.. code-block::
    :number-lines:

    operator (.dot.)
      Dot-product between two vectors a(3) and b(3), v = (a.dot.b).
    note
      - Use parentheses to ensure the correct order of operations.
