#!/bin/bash

#
# This script builds the documentation for my 'qgenrecon' project.  It is set to build
# the documentation for the 'master' branch so that this script may be run as a cron
# job to ensure that the stable version of the code always has up-to-date documentation
# on my server.  
#
# Non-master branches must have their documentation built manually, but this is acceptable
# because non-master branches will be used for development-only.
#
# NOTE: Source code location, along with image location and excluded files/directories,
#	is defined in the doxygen file 'doxyfile.txt'.  This script simply runs DOxygen
#	and copies the result to the appropriate web directory.
#

# Get name of git branch
BRANCHNAME=`git rev-parse --abbrev-ref HEAD`

# Location for html output, using the branch name as a subdirectory
DESTDIR="/srv/www/htdocs/cppdoc/cppclasses/master"


# Checkout 'master' branch
git co master

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

# Remove 'html' directory.  This is to prevent permissions issues if this script is run
# automatically by root.
rm -rf $CURRDIR/"html"

# Checkout branch which was checked-out when this script was created
git co $BRANCHNAME




