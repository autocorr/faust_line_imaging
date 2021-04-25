Developer documentation
=======================
This section is intended to given an overview of the pipeline script
architecture. Contributions are welcome! Simply file a pull request on GitHub
or post an issue for a bug or feature request.

The purpose of the pipeline script is to abstract away the details of imaging in
CASA (calls to ``tclean``, ``immath``, etc.) and replace them with declarations
to higher level pipeline functions specific to FAUST targets and spectral
windows. Interferometric imaging, however, is a time-consuming and often
error-prone process.  Thus in practice the script is meant to be used as a
toolkit that extends the functionality of an interactive session in CASA or a
recipe to be run as a script or batch job. The :doc:`recipes <recipes>` page
details the usage of the :doc:`API <faust_imaging>`. Adding substantial new
functionality, such as new pipeline tasks, will require extending and
understanding the interfaces.

The pipeline code is organized mainly into the categories:

   #. Global variables that contain default parameters.
   #. Utility functions that operate directly on CASA image files. These
      functions are largely independent of project specific information.
   #. Classes that encapsulate information specific to the FAUST project and
      provide methods for the execution of the main pipeline procedures.
      These are composed into the primary user facing class,
      :class:`faust_imaging.ImageConfig`.
   #. Utility functions that perform diagnostic tests, create moment maps,
      and quality assurance plots.

Extending the pipeline mainly requires understanding the four classes:

   #. :class:`faust_imaging.Target`: functionality specific to a FAUST target.
   #. :class:`faust_imaging.Spw`: functionality specific to a FAUST spectral window.
   #. :class:`faust_imaging.DataSet`: functionality specific to an ALMA configuration.
   #. :class:`faust_imaging.ImageConfig`: the composition of the above data classes
      that provides methods to perform the pipeline procedures and ultimately
      call ``tclean``.

Adding new pipeline procedures is thus simply a matter of adding a new method
to implement the functionality. For simple changes or changes intended to be
new default behavior, ``ImageConfig`` should be directly modified with a new
method. For experimental or significant changes that are not intended to be
used as the default, ``ImageConfig`` should be sub-classed and extended,
over-riding the :meth:`faust_imaging.ImageConfig.run_pipeline_tasks` if
necessary. The implementation of pipeline tasks should be factored such that
operations on the "image level" are separated into distinct, small functions
that are composed in the pipeline method/task. If new information or data is
required that is not available in the existing data-classes (``Spw``, etc.),
these should be extended to include it.


