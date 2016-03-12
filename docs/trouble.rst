.. _trouble:

Troubleshooting
===============

Help! ``batman`` is running really slowly - why is this?
--------------------------------------------------------
My first guess is that you're reinitializing the model many times. This is by far the slowest component of ``batman``, because it calculates the optimal step size for the integration starting from a very small value. To speed this up, you can manually set the scale factor for the step size. First check what the optimal step size factor is after you initalize a model with realistic transit parameters:

::
	m = batman.TransitModel(params, t)

	fac = m.fac
	print("stepsize:", fac)

	>>> stepsize: 0.023

Then, you can use set this step size as the default, like so:

::
	m = batman.TransitModel(params, t, fac = fac) 

This will save A LOT of time! The optimal step size will vary slightly depending on the limb darkening and the planet radius, so I would recommend choosing a conservative value for the error tolerance to account for this. 


My light curve has nans in it!
------------------------------
Check and see if your transit parameters are physically realistic. In particular, confirm that the planet's orbit does not pass through the star. (You think this would never happen to you, but don't be too sure! I've been asked about it more than once.)
