.. _tutorial:

Tutorial
============
Here is an example calculation of a light curve with nonlinear limb-darkening (FIXME add other ld models).

First we import ``batman`` and the other modules we'll need (FIXME: make ld a list?):

.. code-block:: python

   import batman
   import numpy as np

Next we initialize the transit model.  We store the transit parameters in a ``TransitParams`` object:

.. code-block:: python

   params = batman.TransitParams()
   params.t0 = 0. 				#transit ephemeris
   params.per = 1.				#orbital period	
   params.rp = 0.1				#planet radius (in stellar radii)
   params.a = 15.				#semi-major axis (in stellar radii)
   params.inc = 1.55				#orbital inclination	
   params.ecc = 0.				#eccentricity	
   params.w = 1.57				#longitude of periastron

and specify the limb-darkening law, the limb darkening coefficients, error tolerance for the model, and the times at which we wish to calculate the model:

.. code-block:: python

   limb_dark = "nonlinear"                 	#limb darkening model
   params.u = np.array([0., 0.7, 0., -0.3]) 	#limb darkening coefficients
   err_max = 0.1                      		#maximum allowed error in light curve (in ppm)
   t = np.linspace(-0.05, 0.05, 1000)    	#times to calculate light curve	

Using these parameters, we initialize the model with the ``TransitModel`` class:

.. code-block:: python

   m = batman.TransitModel(params, t, err_max, limb_dark)

To calculate a model light curve, we use the ``LightCurve`` function: 

.. code-block:: python

   flux = m.LightCurve(params)

Now that the model has been initialized, we can change the transit parameters and calculate a new model like so:

.. code-block:: python
   
   params.rp = 0.11
   new_flux = m.LightCurve(params)

FIXME add link to demo code.


