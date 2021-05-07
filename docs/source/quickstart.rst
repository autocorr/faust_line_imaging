Quick-start guide
=================
While spectral windows (SPWs) are processed in discrete frequency intervals or
"chunks", imaging the full bandwidth cubes can require a large volume of disk
space to store all intermediate products (>100GB for a single SPW and several
TB for a full Setup). After producing the final cubes, the intermediate
products may be removed. Please ensure that the host system has sufficient
resources to run the pipeline. If the host system has limited storage, we
recommend imaging one SPW at a time and removing the intermediate products
after it has finished.


Running the script
------------------
First open CASA from the directory set in ``PROD_DIR`` and then execute the
pipeline script from within CASA:

.. code-block:: python

   # using an absolute path
   execfile('/<PATH>/<TO>/faust_imaging.py')
   # or using a relative path
   execfile('../casa_scripts/faust_imaging.py')

This ``execfile`` command can also be performed in scripts that are themselves
executed with ``execfile`` in CASA.


Running the pipeline
--------------------
To run the imaging pipeline for one target and SPW, create an instance of the
class :class:`faust_imaging.ImageConfig` (these links point to the API
documentation which further describe the calling convention) and call the
:meth:`faust_imaging.ImageConfig.run_pipeline` method. Images are referenced by
their field name and a unique label for the SPW based on the molecular tracer
and rest frequency. To see these names, one can print the global variables:

.. code-block:: python

   # See all SPW labels
   print(ALL_SPW_LABELS)

   # See all field names
   print(ALL_FIELD_NAMES)

Based on the above, to image the CS (5-4) transition in Setup 2 for source
CB68, one would run:

.. code-block:: python

   config = ImageConfig.from_name('CB68', '244.936GHz_CS')
   config.run_pipeline()

The default parameters will use a Briggs robust *uv*-weighting of 0.5 and
jointly deconvolve all array configurations (12m & 7m). We recommend to image
one line first to test that the pipeline works and produces sensible results
before moving batched processing. The final pipeline products will have a
suffix "_clean" before the file extension (e.g., ".image"). In the above
example, the final image (primary beam corrected) will be written to:

.. code-block:: bash

   # Path relative to where CASA is run (also the value of `PROD_DIR`)
   images/CB68/CB68_244.936GHz_CS_joint_0.5_clean.image.pbcor.common.fits

   # Files are named according to the convention:
   #   <FIELD>_<SPW_LABEL>_<ARRAY>_<WEIGHTING>_<SUFFIX>.<EXT>
   # where
   #   FIELD     : source/field name
   #   SPW_LABEL : rest-frequency and molecular tracer based name
   #   ARRAY     : array configurations used; can be 'joint', '12m', '7m'
   #   WEIGHTING : uv-weighting applied, e.g. 'natural' for natural weighting
   #               or '0.5' Briggs robust of 0.5.
   #   SUFFIX    : name to distinguish different files, 'clean' is used for
   #               the default final pipeline products. Intermediate products
   #               will also exist with suffixes including 'dirty', 'nomask',
   #               etc.
   #   EXT       : CASA image extension name, e.g. '.image' or '.mask'

Please see the :doc:`API documentation <faust_imaging>` or docstring for
further configuration information. For examples of more advanced uses
of the pipeline please refer to the :doc:`Cookbook <recipes>`.


Running the serial pipeline
---------------------------
To run the imaging pipeline for all SPWs of a target one at a time, use the
:func:`faust_imaging.run_pipeline` function. Note that this is likely to be
slow for most systems and could take approximately one week or longer to run.

.. code-block:: python

   # for just Setup 2
   run_pipeline('CB68', setup=2)
   # for all setups (default)
   run_pipeline('CB68', setup=None)


Running the parallel pipeline
-----------------------------
To run the pipeline in parallel, please refer to the :any:`ParallelCasa`
section in the Cookbook. Example scripts are included for imaging a single
SPW in parallel and also imaging all of the SPWs for a setup in parallel.
On the NRAO NM postprocessing cluster, typical run-times are a few hours
when imaging a single SPW in parallel and a few days for imaging all SPWs
of a setup.


.. _Quick Test:

A quick test
------------
To quickly test that the pipeline works from (approximately) beginning to end,
a single small frequency interval can be imaged by itself:

.. code-block:: python

   # Use a field and SPW with measurement-set data available on your host
   full_config = ImageConfig.from_name('CB68', '244.936GHz_CS')
   chunked_configs = full_config.duplicate_into_chunks(nchunks=100)
   first_chunk = chunked_configs[0]
   first_chunk.run_pipeline()

Further details are provided in the :doc:`Cookbook <recipes>`.


