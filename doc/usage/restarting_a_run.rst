Restarting a run
================

It is possible to restart a run that has terminated abrutly, or that
has terminated succesfully, but needs to be further extended. This is
done by rerunning the job while providing the path of the lastly saved
binary state file (see :doc:`/usage/output_files`) as the parameter to
the option ``--restart`` (see :doc:`/usage/command-line`).

When restarting a run, the input file can be altered so as to override
the original input parameters stored in the state file. This can be
used, for example, to allow for more grid-shifting iterations to be
executed, after the maximum number has already been reached. The
parameters of which the values in the input file override those in the
state file are indicated as *mergeable* in their corresponding
descriptions in the section :doc:`/usage/input_file`. For a quick
reference, they are also listed below:

.. include:: blocks/mergeables.rst
