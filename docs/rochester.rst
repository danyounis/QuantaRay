Documentation
=============
Module: Rochester Potential ``(rochester.f08)``
-----------------------------------------------
Revised 9/14/2025

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
