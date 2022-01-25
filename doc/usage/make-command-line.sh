#!/bin/bash

cat > command-line.rst << EOF
############
Command-line
############

.. code-block::

$(gmak -h | sed "s/^/    /")
EOF
