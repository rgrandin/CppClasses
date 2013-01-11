#!/bin/bash

#
# This script builds the documentation for my 'qgenrecon' project.  The name of the active
# branch within the git repository is obtained and used as a subdirectory within the web
# folder containing documentation for 'qgenrecon'.
#
# NOTE: Source code location, along with image location and excluded files/directories,
#	is defined in the doxygen file 'doxyfile.txt'.  This script simply runs DOxygen
#	and copies the result to the appropriate web directory.
#

# Get name of git branch
BRANCHNAME=`git rev-parse --abbrev-ref HEAD`

# Location for html output, using the branch name as a subdirectory
DESTDIR="/srv/www/htdocs/cppdoc/cppclasses/"$BRANCHNAME

# Get current directory.
CURRDIR=`pwd`

# Delete 'html' directory to ensure that the built documentation is current.
rm -rf $CURRDIR/"html"

# Run DOxygen on copied files.  Documentation is generated in local directory
# and then copied to desired location.
doxygen doxyfile.txt

# Make the destination directory, if it doesn't already exist
mkdir -p $DESTDIR

# Copy 'html' output to destination directory.
cp -r $CURRDIR/"html" $DESTDIR





