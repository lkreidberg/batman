.. batman documentation master file, created by
   sphinx-quickstart on Tue Jun 16 12:35:22 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. |br| raw:: html

   <br />

batman: |br| Bad-Ass Transit Model cAlculatioN
=========================================

The ``batman`` package for Python makes super fast calculation of transit light curves easy.  

The package supports quadratic and nonlinear limb-darkening, as well as custom user-specified profiles. The lion’s share of the computation is done with C extension modules, so it is very fast. ``batman`` is parallelized with OpenMP and compatible with Python 2-3.

Contents:

.. toctree::
   :maxdepth: 2

   installation
   quickstart
   tutorial
   api
   acknowledgements

Release Notes
-------------
.. include::  ../changelog.rst


