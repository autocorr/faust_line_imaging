User Guide
==========
Welcome to the pipeline User Guide! This guide details how to (1) run the
pipeline interactively for combinations of targets and spectral windows (SPWs),
run batch scripts for parallel processing, and how to trouble-shoot the
pipeline products and results. A table of contents tree may be found in the
left side-bar.


Note on system resources
------------------------
While spectral windows are processed in discrete frequency intervals or
"chunks", imaging the full bandwidth cubes can require a large volume of disk
space in order to store all intermediate products. For example, the products
can be greater than 100GB for a single SPW and several TB for a full Setup.
After producing the final cubes, the intermediate products may be removed.
Please ensure that the host system has sufficient resources to run the
pipeline. If the host system has limited storage, we recommend imaging one SPW
at a time and removing the intermediate products after it has finished.


Importing the pipeline
----------------------
After ensuring that the paths are configured correctly (i.e., ``DATA_DIR`` points
to the directory containing the measurement sets and CASA is started from the
directory set in ``PROD_DIR``). Start CASA interactively from a shell prompt
and execute the pipeline script using ``execfile`` (see the :doc:`Installation
<environment>` page). For an example ``PROD_DIR="/mnt/scratch/faust"`` and
the pipeline script located in ``/mnt/scratch/faust/faust_line_imaging/faust_imaging.py``
then an interactive session can be started by running:

.. code-block:: bash

   $ cd /mnt/scratch/faust  # for example
   $ casa

and then from the interactive CASA prompt:

.. code-block:: python

   CASA <1>: execfile('faust_line_imaging/faust_imaging.py')


Names and labels
----------------
Targets are referred to by their field names, e.g. "CB68" or "BHB07-11", as
they appear in the measurement sets. All valid field names can be read from
the global variable ``ALL_FIELD_NAMES``:

.. code-block:: python

   print(ALL_FIELD_NAMES)

SPWs are assigned unique labels based on the principle molecular tracer
targeted in the window. These labels are prefixed by the truncated rest
frequency of the transition and followed by the simple for the molecule. The CS
(5-4) transition in Setup 2 has the label "244.936GHz_CS". All valid SPW labels
can be read from the global variable ``ALL_SPW_LABELS``:

.. code-block:: python

   print(ALL_SPW_LABELS)

To image narrow windows around secondary lines without imaging the full bandpass
of a SPW, see the Cookbook section :any:`Imaging cut-out velocity windows`.


Running the pipeline
--------------------
The primary user interface for running the pipeline tasks is through the Python
class :class:`faust_imaging.ImageConfig`, specifically the constructor
:meth:`faust_imaging.ImageConfig.from_name`. These and similar links point to
the API documentation that further describe the calling convention and
arguements functions and classes take. Most of this information can also be
found in the docstrings as well, which can be read by running ``help <ITEM>``
or ``<ITEM>?`` from the CASA prompt. Once an instance of ``ImageConfig`` is
created and set with the desired options, the pipeline can be be run using the
:meth:`faust_imaging.ImageConfig.run_pipeline` bound method. As a minimal
example:

.. code-block:: python

   config = ImageConfig.from_name('CB68', '244.936GHz_CS')
   config.run_pipeline(ext='clean')
   # The final products will be given a suffix based on `ext` above,
   # the default is "clean". Different names may be used to avoid over-
   # writing existing files.

The default parameters will use a Briggs robust *uv*-weighting of 0.5 and
jointly deconvolve all array configurations (12m & 7m). We recommend imaging
one line first to test that the pipeline works and produces sensible results
before moving to batched processing. The final pipeline products will have a
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

   # Other non-standard image extensions produced will include '.common' for
   # images that have been smoothed to a common beam resolution, '.hanning' for
   # hanning smoothed images, and '.pbcor' for images corrected for attenuation
   # by the primary beam.

Please see the :doc:`API documentation <faust_imaging>` or docstring for
further configuration information. For examples of more advanced uses
of the pipeline please refer to the :doc:`Cookbook <recipes>`.


Quality assurance plots and moment maps
---------------------------------------
After the pipeline has been run, the next steps are to validate the results by
creating quality assurance plots for visual inspection and moment maps. The QA
plots generate channel maps of the restored image and residual with the
clean-mask overplotted. Only channels containing significant emission are plotted
(regardless of whether the emission is masked). The default threshold to show
such channels is 6-times the full-cube RMS.  To create the quality assurance
plots call the function :func:`faust_imaging.make_all_qa_plots` for the desired
field and extension (e.g., "clean" as used above).

.. code-block:: python

   make_all_qa_plots('CB68', ext='clean', overwrite=False)

The ``overwrite=False`` keyword argument ensures that QA plots are only
generated for images that do not already exist, so this function can be safely
called after new pipeline jobs have been run. Plots will be written to
``plots/`` or the value of ``PLOT_DIR``. Note that the creation of the plots is
implemented inefficiently with ``matplotlib`` and ``imshow``, and thus creating
the plots for the Setup 3 SPWs may require >50 GB of memory.

.. image:: images/cb68_qa_image_238.png
   :width: 300

.. image:: images/cb68_qa_residual_238.png
   :width: 300

The above figures show channel index number 238 of CB68 CS (5-4) for the
restored image (**left**) and the residual image (**right**). For the restored
image, color-scale ranges from -3 to 10 times the full-cube RMS and the filled
contours are shown at increments of 10, 20, 40, and 80 times the RMS.  The cyan
contour shows the clean mask. For the residual image the color scale ranges
from -5 to 5 times the RMS (negatives shown in blue) and the black contour
shows the clean mask. The tick-marks show increments of 5 arcsec, the dashed
line shows the half-power beamwidth of the primary beam, and imaged out to the
20%-power point of the primary beam.

Further details on the QA plots may be found in the Cookbook section
:any:`QA Plots`.

Moment maps may be generated by running the :func:`make_all_moment_maps`
function for the desired field and extension (e.g., "clean" as used above).

.. code-block:: python

   make_all_moment_maps('CB68', ext='clean')

Maps based on the integrated intensity ("mom0"), maximum or peak intensity
("max"), centroid velocity ("mom1"), and velocity dispersion ("mom2") will be
written to the ``moments/`` directory or the value of ``MOMA_DIR``. The images
can be inspected with your FITS viewer of choice. The figure below shows an
example matplotlib visualization of the moments for CB68 CS (5-4).

.. image:: images/cb68_moments.png
   :width: 640

The above figures can be generated using the ``util/moment_plotting.py`` script
under Python **v3** (currently undocumented; requires packages numpy, scipy,
skimage, matplotlib, aplpy, radio_beam, and astropy).


Trouble-shooting
----------------
Having made the deconvolved image products and the quality assurance plots, the
next step is to inspect the results and resolve whether they are satisfactory
for the science-goals of the Source Team.  A few common cases where the
deconvolution and/or masking have produced undesirable results are detailed
below.

Extended negative emission
~~~~~~~~~~~~~~~~~~~~~~~~~~
If the emission is extended and has negative-intensity artifacts or "bowls" due
to missing short spacings, has ``tclean`` masked any of these negative
features? It is undesirable to include these artifacts in the source model.
If they are included, the auto-masking parameters may be tuned to limit
the masking of negative emission.

.. code-block:: python

   config = ImageConfig.from_name('CB68', '244.936GHz_CS')
   config.autom_kwargs['negativethreshold'] = 8  # the default is 7
   config.run_pipeline()

True absorption does frequently occur, however, towards the bright and compact
continuum emission the central protostellar source(s). Because the visibility
data is continuum subtracted, this absorption will appear negative in the
restored images.  This absorption should be masked and cleaned.

Overly-aggressive clean masks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Does the automated clean-mask appear to be overly aggressive and include large
areas of what does not appear to be genuine emission?  This effect has been
known to appear in earlier iterations of the pipeline for certain fields.  The
nature of the auto-multithresh algorithm gives these spurious masks an "amoeba"
or "algae" like appearance, as can be seen in the following figure:

.. image:: images/amoeba_example.png
   :width: 400

In some circumstances, all pixels in a channel may be included in the mask.
Note that such cases will appear to have no mask when using the ``casaviewer``
to plot a contour-diagram.  This effect seems to be largely mitigated with the
latest set of default parameters, but careful attention should be paid in case
it appears.

Spurious and overly-aggressive masking may change the noise statistics and
include artifacts in the source model and should be re-imaged.  The most
straightforward solution is to raise the significance threshold used to "grow"
the mask.

.. code-block:: python

   config = ImageConfig.from_name('CB68', '244.936GHz_CS')
   config.autom_kwargs['lownoisethreshold'] = 2.0  # the default is 1.5
   config.run_pipeline()

Significant uncleaned emission
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If the automated masking appears to have left significant emission unmasked
and thus uncleaned. This can frequently be diagnosed in the QA plots of the
residual image. The investigator may use their discretion to decide whether
such emission produces adverse affects and should be cleaned.  Multiple methods
exist to fix such images without performing the full pipeline over again.
Namely, the final clean may be restarted with:

   #. Using auto-multithresh but with a lower 'lownoisethreshold'.
   #. Using auto-multithresh and manually adding regions to the existing mask.
   #. Without using auto-multithresh and manually adding regions to the
      existing mask.

The pipeline processes discrete image "chunks" in frequency to improve
performance and ease memory constraints. Restarting thus requires operating
on the chunk containing the offending emission. In the following example,
channel index number 238 is insufficiently cleaned and the offending chunk
is restarted with the interactive cleaning.

.. code-block:: python

   full_config = ImageConfig.from_name('CB68', '244.936GHz_CS')
   chunked_configs = full_config.duplicate_into_chunks()
   problematic_config = chunked_configs.get_chunk_from_channel(238)
   problematic_config.clean_line(ext='clean', restart=True, interactive=True)
   # ^ The casaviewer will appear for manual masking. Identify the channel
   #   with the offending emission (the channel indices will now be of the chunk)
   #   and draw an addition to the mask. Often times it suffices to select
   #   the "blue rightward arrow" icon immediately if the emission is faint.
   chunked_configs.postprocess(ext='clean')

More information on manually restarting one chunk is described in the Cookbook
:any:`Restarting one chunk` section.

Inconsistent masking from varying noise
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The sensitivity as a function of frequency for some SPWs is affected by
atmospheric or telluric lines. Examples include the "231.221_13CS" and
"231.322_N2Dp" SPWs in Setup 1. An atmospheric ozone feature between these two
windows increases the RMS by about 20% towards the respective band edge.  In
some circumstances, the use of a single RMS can lead to over-masking of many
small noise spikes near the band edge. If this is the case, then using smaller
image-chunk sizes should give more uniform results.

Divergences or negative edge-features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
It is a known issue that the multiscale clean implementation in CASA can
introduce instability when using clean masks. In some circumstances
``tclean`` can diverge at the edge of the clean mask or primary beam mask
and insert spurious positive-intensity features into the model. These
features are usually on large scales (often similar to the ACA synthesized
beam) and produce strong negative-intensity features in the restored
image.

The default parameters have been found to largely stabilize ``tclean`` by
slowing the rate of convergence in the minor cycle. If these divergences
appear, try running the pipeline with a lower ``gain`` and higher
``cyclefactor``:

.. code-block:: python

   config = ImageConfig.from_name('CB68', '244.936GHz_CS')
   config.gain = 0.03  # default 0.05
   config.cyclefactor = 2.5  # default 2.0
   config.run_pipeline()


Running the parallel pipeline
-----------------------------
To run the pipeline in parallel, please refer to the :any:`Parallel CASA`
section in the Cookbook. Example scripts are included for imaging a single
SPW in parallel and also imaging all of the SPWs for a setup in parallel.
On the NRAO NM postprocessing cluster, typical run-times are a few hours
when imaging a single SPW in parallel and a few days for imaging all SPWs
of a setup.


