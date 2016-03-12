.. batman documentation master file, created by
   sphinx-quickstart on Tue Jun 16 12:35:22 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. |br| raw:: html

   <br />

batman: |br| Bad-Ass Transit Model cAlculatioN
==============================================

Welcome to the documentation for ``batman``, a Python package for fast calculation of exoplanet transit light curves.  The package supports calculation of light curves for any radially symmetric stellar limb darkening law, using a new integration algorithm for models that cannot be quickly calculated analytically. 

In typical use, ``batman`` can calculate a million model light curves in well under 10 minutes for any limb darkening profile.

A `paper <http://arxiv.org/abs/1507.08285>`_ describing ``batman`` is up on arXiv; please cite it if you make use of the package!

Contents:

.. toctree::
   :maxdepth: 2

   installation
   quickstart
   tutorial
   api
   trouble 
   acknowledgements

Release Notes
-------------
.. include::  ../changelog.rst


