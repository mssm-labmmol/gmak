##########
Input File
##########

Except for restart jobs, ``gmak`` initializes its internal objects by
reading data from an input text file. The lines in this file are
organized into blocks. Each block begins with a line ``$blockname``,
where ``blockname`` is the name of the block, and ends with a line
``$end``.  Most input parameters are set within blocks that identify
their context in the program.  Those that are not set inside any
blocks are referred to as global input parameters.

Parameters are set using lines with syntax::

    PARAMETER VALUE...

where ``PARAMETER`` is the name of the parameter, ``VALUE`` is a string
with no white spaces (referred to as a value token), and the ellipsis
``...`` indicates that there may be a list of such tokens in the line.
Tokens can be separated from the parameter name and from other tokens by
any number of spaces.

The tokens can be categorized into the following types:

a. Numerical

b. File/directory path

c. Object reference (the name of an object previously initialized in
   the input file)

d. Mathematical expression

e. String (when none of the above apply)

Any line that starts with the ``#`` character is considered a comment.

Global input parameters
=======================

Global input parameters should appear before any blocks.  They either
initialize data or set behaviors that belong to no particular block
context.

The table below lists the known global input parameters, as well as
the expected type of value tokens.

.. include:: blocks/toc.rst

Blocks
======

Blocks set a contextual scope for initializing data and behaviors.
Typically, each instance of a block translates as the initialization
of an object of which the type is tightly associated with the block
name. For example, a ``$protocol`` block sets up an instance of a 
class that implements a :doc:`protocol </overview/protocols>`.

Blocks may depend on each other via object-reference tokens or due to
the logic of interpretation of the input file. In the latter case, the
dependency is classified as *intrinsic*, and the dependent block must
appear after all instances of the one on which it depends. In the
former case, the block that defines the object being referred to
must precede the block that refers to it. 

The table below lists the known blocks and, if existing, their
intrinsic dependencies on other blocks. 

.. include:: blocks/blocks.rst

.. include:: blocks/toctree.rst

