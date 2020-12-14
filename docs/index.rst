.. FEA documentation master file, created by
   sphinx-quickstart on Tue Dec  5 12:39:10 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Flux Envelope Analysis
===============================

This is a Python package for calculating a reduced-dimensional solution space from a larger linear program using optimization. It was designed for solving 2D and 3D solution spaces from large linear programs as utilized in constraint-based genome-scale metabolic modeling.

It leverages `optlang <https://github.com/opencobra/optlang>`_ to interface with linear program solvers. This is also the backend that is utilized by `cameo <http://cameo.bio/>`_.

It includes some :class:`matplotlib` plotting functions for generating basic 2D and 3D graphs from the reduced solution spaces (and even 4D if you're willing to use color to represent the fourth dimension).

For more details on the implementation, read the `full manuscript <../_static/FEA_Manuscript.pdf>`_. For examples and current issues, see the `GitHub page <https://github.com/long-m-r/fea>`_.

Documentation
=======================

The main entry point into FEA is the flux_envelope_analysis function documented here. Other documentation is linked below.

.. automodule:: fea
   :members:

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   plotting
   internal classes

Index
=======

* :ref:`genindex`
* :ref:`modindex`