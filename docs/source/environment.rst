Pipeline installation and setup
===============================
The pipeline comprises a single Python file to be run under CASA v5. The
pipeline script does not require installation per se, but must be executed when
starting a new CASA process or session. The instructions below detail how to
download the script, setup the correct pathes and directories, and execute the
script.


CASA compatibility
-----------------
The pipeline is to be run using CASA *v5* (e.g., v5.6, v5.7), which will be the
last version of CASA to use Python 2 and the included CASA Python distribution.
CASA v6 and Python v3 are not supported at this time, but may be in the future.


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
the measurement sets and write out the images. The following global paths
are to be set by the user:

.. code-block:: python

   DATA_DIR  # where the measurement sets are stored
   PROD_DIR  # where CASA will be run and the output products stored

Edit these paths in ``faust_imaging.py`` for the correct values on your host
system. Directories in ``PROD_DIR`` will be automatically created to store
images, moment maps, and plots.

Visualizing the directory structure as a tree would have files organized as
such

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

Some CASA tasks do not work with absolute paths, so please note that the
pipeline script *must* be run from the directory specified by ``PROD_DIR``.
The measurement set files should be present in the directories:

.. code-block:: bash

   $DATA_DIR/<TARGET>-Setup1/
   $DATA_DIR/<TARGET>-Setup2/
   $DATA_DIR/<TARGET>-Setup2/

where ``<TARGET>`` is the FAUST target field name, e.g. "CB68" or "L1527".
The calibrated measurement sets may be downloaded from RIKEN. The names may be
found in the ``ALL_TARGETS`` dictionary, I retrieved these values from the
proposal, so they may be inconsistent for some targets. The paths above may
modified directly by editing the format string attribute ``DataSet.ms_fmt``.


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


Next steps
----------
Congratulations! Now that the environment is setup, please now refer to the
:doc:`User Guide <userguide>` or click the "Next" button for instructions on
running the pipeline.


