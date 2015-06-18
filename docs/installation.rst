.. _installation:

Installation
============

Download the source directory from https://github.com/lkreidberg/batman/archive/master.zip.

Unpack the distribution with:

::

   $ unzip batman-master.zip

To install, navigate to the source root directory and run the setup script:

::

   $ cd batman-master
   $ sudo python setup.py install

Finish by cleaning up:

::
   
   $ cd ..
   $ sudo rm -rf batman-master

Now ``batman`` is installed and ready to use! To verify that the installation is working properly, run a few tests with:

::

   $ python -c 'import batman; batman.test()'


