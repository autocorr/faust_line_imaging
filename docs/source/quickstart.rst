Quick-start guide
=================
First note that many of these image products are very large take a great deal
of memory to run (> 60GB). Further, when imaging full bandwidth cubes, the
resulting products can require a large volume of disk space (potentially >
100GB for a single SPW). After producing the final cubes, many of these
products can be removed however. Please ensure that the host system has
sufficient resources.

First open CASA from the directory set in ``PROD_DIR`` and then execute the
pipeline script from within CASA:

.. code-block:: python

   # using an absolute path
   execfile('<PATH>/<TO>/faust_imaging.py')
   # or using a relative path
   execfile('../casa_scripts/faust_imaging.py')

To run the imaging pipeline for all SPWs of a target, use the ``run_pipeline``
function now in the global scope.

.. code-block:: python

   # for just Setup 1
   run_pipeline('CB68', setup=1)
   # for all setups (default)
   run_pipeline('CB68', setup=None)

A list of uv-weightings to use and whether to image the fullcube or a windowed
sub-cube may be passed.

.. code-block:: python

   # image both Briggs robust=0.5 and Natural uv-weightings
   # image a 20km/s window centered on the priamry line instead of all channels
   run_pipeline('CB68', setup=1, weightings=(0.5, 'natural'), fullcube=False)

Please see the :doc:`API documentation <faust_imaging>` or docstring for
further configuration information. For imaging individual SPWs or more
fine-grained control of the pipeline, please see the following :doc:`cookbook
<recipes>` section.
