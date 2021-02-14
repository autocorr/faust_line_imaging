Pipeline Cookbook
=================
The heuristics adopted in the pipeline appear to work well in the majority of
cases, but there are invariably cases not adequately treated by the defaults
and thus require custom processing. To aid in the processing of a FAUST target,
the pipeline provides helper classes that abstract away the book-keeping aspects
into individual tasks to be configured depending on the requirements of a
particular field and SPW.

Note that because the pipeline is computationally intensive, even imaging
a single SPW can take a day or more when CASA is run using single threaded
execution, and still a substantial amount of time when executed in parallel.

Running the default pipeline
----------------------------
The :class:`faust_imaging.ImageConfig` class provides the primary interface to
``tclean`` within CASA and encapsulates properties specific to a field, SPW,
array configuration, and desired ``tclean`` parameters. Please refer to the
:doc:`API Documentation <faust_imaging>` and the docstring for additional
information on the calling convention of this class.

To run all tasks of the pipeline with default parameters, first create an
instance of the :class:`faust_imaging.ImageConfig` class and use the
:meth:`faust_imaging.ImageConfig.run_pipeline` method.

.. code-block:: python

   config = ImageConfig.from_name('CB68', '244.936GHz_CS', weighting=0.5)
   config.run_pipeline()

The full list of SPW labels may be found in the `ALL_SPW_LABELS` variable.
The above command will generate the default pipeline image products for
the target field CB68, for the CS (5-4) line of Setup 2, using a Briggs
robust factor of ``0.5``. The full calling convention is:

.. code-block:: python

   ImageConfig.from_name(
           '<FIELD>',  # field name, e.g., "CB68". See the global var `ALL_FIELDS`
           '<LABEL>',  # SPW label, e.g., "244.936GHz_CS"
           weighting=0.5,
           # ^ Use a number for Briggs robust or a string recognized by tclean,
           # e.g., 'uniform' or 'natural'. The default is 0.5.
           fullcube=True,
           # ^ Image the full bandwidth of the SPW or a narrow 20km/s window
           # centered on the primary line of interest. The default is True.
           kind='joint',
           # ^ What combination of array configurations to use. Possible values
           # include ('joint', '12m', '7m'). The default is 'joint'.
   )

Pipeline task description
-------------------------
The :meth:`faust_imaging.ImageConfig.run_pipeline` method discussed above is
primarily a wrapper for calling all of the pipeline tasks in sequence using the
default parameters. For custom imaging, it is recommended that users create a
recipe using the underlying tasks.  As an example, this is the code that is
executed by default:

.. code-block:: python

   config = ImageConfig.from_name(...)
   ext = 'clean'  # extensionn name for final cleaned products
   config.make_dirty_cube()
   config.clean_line_nomask(sigma=4.5)
   config.make_seed_mask(sigma=5.0)
   config.clean_line(mask_method='seed+multithresh', ext=ext)
   config.postprocess(ext=ext, make_fits=True)

``TODO`` describe each task and the calling convention detail.

Restarting ``tclean`` and manual masking
----------------------------------------
``TODO``

Parallel image processing
-------------------------
The computational efficiency for imaging a single SPW can be improved by
a factor of approximately two using CASA run with MPI (i.e., ``mpicasa``).
The pipeline script will automatically recognize from the environment whether
it is being run with ``mpicasa`` and use ``parallel=True`` accordingly in
``tclean``. At the current time the preferred pipeline masking method of
``mask_method='seed+multithresh'`` is not currently supported in parallel.
Imaging may be run in parallel when using ``mask_method='auto-multithresh'``,
although the former method is preferred.

The majority of the computational run-time cost of the imaging derives from
the conservative ``cyclefactor`` and ``gain`` used with ``tclean``. These
parameters are selected to avoid divergences that are common when cleaning
extended emission with large scales using multiscale clean and clean
boxing/masking. If the emission is relatively compact in the field and
SPW of interest, these parameters can be relaxed to increase performance:

.. code-block:: python

   config.gain = 0.1         # default 0.05
   config.cyclefactor = 1.0  # default 2.0
   config.run_pipeline()

