Bibliography
============

meq is an implementation of published methods, and this documentation cites
them where they are used rather than collecting them into a reading list. This
page is the resolved bibliography.

Two of these are not background reading — they *are* the method.
:cite:t:`SanchezVizuetSolano2019` fixes the spaces, the numerical flux and the
stabilisation; :cite:t:`SanchezVizuet2020adaptive` adds adaptivity, the error
estimator and the treatment of the curved boundary, and supplies meq's headline
benchmark. Between them they determine nearly every choice in ``src/meq``. Read
the first completely and work from the second thereafter.

``refs/Refs.md`` in the source tree is the same list with the DOI to fetch each
PDF by, plus annotations recording what each paper contributes and — for the
ones that have been read against meq's own formulation — what it gets wrong.
The PDFs themselves are not tracked, being publisher material.

.. note::

   **Three of the papers meq depends on carry errata that this project found by
   checking them**, and all three are of the kind that converge beautifully to
   the wrong answer: a sign on a source term, a duplicated coefficient in a
   published table, and two reversed signs in a derivative that neither of that
   paper's own benchmarks can detect. Where meq's documentation departs from a
   published equation, it says so and says why. The individual cases are in
   :ref:`formulation-soloviev-sign`, :ref:`testing-fixtures` and
   :ref:`rotation-errata`.

   This is not a criticism of the papers. It is the reason every fixture in
   ``tests/analytic/`` recomputes :math:`\dstar\psi` by finite differences and
   asserts it against the source it hands the solver, instead of trusting a
   transcription.

.. bibliography::
   :all:
