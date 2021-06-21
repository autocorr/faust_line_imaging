Pipeline installation and setup
===============================
The pipeline comprises a single Python file to be run using `CASA <https://casa.nrao.edu/>`_.
The pipeline script does not require installation per se, but must be executed
when starting a new CASA process or session. The instructions below detail how
to download the script, setup the correct pathes and directories, and execute
the script.


CASA compatibility
-----------------
The pipeline script can be used with the "monolithic" versions of CASA v5 (Python 2) or
CASA v6 (Python 3). These are the versions of CASA downloaded from the main website.
The pipeline has been primarily tested under CASA v5.6 (the current Calibration
Pipeline release), but has been tested to at least work under CASA v6.2.
The "Standard Products" to be distributed to the ALMA Archive are made using CASA v5.6.

The script may also be run using the
`casatasks <https://open-bitbucket.nrao.edu/projects/CASA/repos/casa6/browse/casatasks>`_
and `casatools <https://open-bitbucket.nrao.edu/projects/CASA/repos/casa6/browse/casatools>`_
Python 3 library modules.
The latter may be installed into the user's Python environment by running
``pip install --user casatasks casatools``,
preferably within a project-specific virtual environment.


Downloading the pipeline
------------------------
The pipeline comprises a single Python source file, ``faust_imaging.py``, that
may be downloaded either by downloading the repository from `GitHub
<https://github.com/autocorr/faust_line_imaging>`_ as a Zip archive or by
cloning the repository using Git:

.. code-block:: bash

   git clone https://github.com/autocorr/faust_line_imaging
   # this will create the directory "faust_line_imaging" in the current
   # working directory.

Updates to the script from the main repository can be automatically merged by
running ``git pull`` from the repository. To modify the source and have those
changes reflected in the main repositry, please `fork
<https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/cloning-and-forking-repositories-from-github-desktop#forking-a-repository>`_
the above repository on
GitHub and file a `pull request
<https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request-from-a-fork>`_.


Paths and directory setup
-------------------------
The pipeline script currently requires an explicit directory structure to read
the measurement sets and write out the images. Essential global paths are set
by the user in the configuration file ``faust_imaging.cfg``. The contents of
the template file included in the source repository are written below:

.. code-block:: cfg

   # This is the user configuration file for the FAUST archival imaging pipeline.
   # Edit this file with settings appropriate for your host system and place it in
   # the directory specified by PROD_DIR.

   [Paths]
   # Specify where the measurement set files/directories are located. The MSs
   # for a field should follow the directory structure:
   #   <DATA_DIR>/<FIELD>-Setup<N>/*.ms
   # e.g.,
   #   /scratch/faust/CB68-Setup1/*.ms
   DATA_DIR=/lustre/aoc/users/USER/FAUST/2018.1.01205.L/completed_SBs/

   # Specify where CASA is to be run and also the base-path where products
   # are to be written.
   PROD_DIR=/lustre/aoc/users/USER/faust/faust_alma/data/

Copy the template ``faust_imaging.cfg`` file from the repository into
the directory where CASA is to be run from (``PROD_DIR``) and edit
the paths appropriately for your host system. This file will be read
when the pipeline script is executed in CASA. Sub-directories for images
cubes, moments, and plots will be automatically generated.

Visualizing the directory structure as a tree, the file hierarchy would
be organized:

.. code-block:: none

   /<PATH>/<TO>/DATA_DIR/   <- location on user's system
     - CB68-Setup1          <- MS files per target per setup
     - CB68-Setup2
     - ...
   /<PATH>/<TO>/PROD_DIR/   <- location on user's system
     - images/              <- these will be made automatically
       + CB68/              <- image cubes and products per target
       + L1527/
       + ...
     - moments/
       + CB68/              <- moment maps per target
       + L1527/
       + ...
     - plots/               <- diagnostic and QA plots

Some CASA tasks do not work using absolute paths, so please ensure that the
pipeline script is run from the directory specified by ``PROD_DIR``.  The
measurement set files should be present in the directories:

.. code-block:: bash

   $DATA_DIR/<TARGET>-Setup1/
   $DATA_DIR/<TARGET>-Setup2/
   $DATA_DIR/<TARGET>-Setup2/

where ``<TARGET>`` is the FAUST target field name, e.g. "CB68" or "L1527".
The calibrated measurement sets may be downloaded from RIKEN. The names may be
found in the ``ALL_TARGET_NAMES`` global variable, I retrieved these values
from the proposal, so they may be inconsistent for some targets. The paths
above may modified directly by editing the format string attribute
``DataSet.ms_fmt``.


Executing the script
--------------------
Ensure that the script can be properly executed from within CASA by starting
CASA from the directory set in ``PROD_DIR``. The pipeline script can then be
executed and all functions/symbols brought into scope using the ``execfile``
command:

.. code-block:: python

   # using a relative path
   execfile('../<Path>/faust_line_imaging/faust_imaging.py')
   # or alternatively using an absolute path
   execfile('/<PATH>/<TO>/faust_line_imaging/faust_imaging.py')

This ``execfile`` command needs to be run whenever starting CASA or when
the pipeline script source code is modified.  Note that the ``execfile``
command can also be performed within "recipe" scripts that are themselves run
with ``execfile`` in CASA.


Notes for MacOS users
---------------------
MacOS users will likely need to run the command ``ulimit -Sn 8000`` from the
shell before starting CASA. This command increases the maximum number of files
that may be opened at once.


Next steps
----------
Congratulations! Now that the environment is setup, please now refer to the
:doc:`User Guide <userguide>` or click the "Next" button for instructions on
running the pipeline.


