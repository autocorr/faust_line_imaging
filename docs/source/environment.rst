Filesystem and CASA environment
===============================
The pipeline script ``faust_imaging.py`` is to be run using CASA v5 (e.g.,
v5.6), which will be the last version of CASA to use Python 2 and the included
CASA Python distribution. The pipeline may support CASA v6 and Python v3 in the
future, but does not currently at this time.

The pipeline script currently requires an explicit directory structure to read
the measurement sets and correctly write out the images. Note the following
global paths to be set by the user in particular:

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

Now that the environment is setup, please see the :doc:`quick-start guide
<quickstart>` or click the "Next" button for running the pipeline.
