#!/bin/bash

#
# This script accepts a single file as an input argument and writes the git
# repository hash for that file to the terminal.
#

HASH=`git hash-object $1`

echo $HASH
