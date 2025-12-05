#!/bin/bash
# A script to compile the AOMC Fortran project.

# Exit immediately if a command exits with a non-zero status.
set -e

# Compiler
FC=gfortran

# List of source files, now in the src/ directory
MODULE_FILES="src/global.f90 src/modules.f90"
SOURCE_FILES="src/geom2.f90 src/interface.f90 src/light_internal.f90 src/logbin.f90 src/mc.f90 src/vsf.f90 src/water.f90"
ALL_OBJECTS=""

echo "Compiling modules..."
$FC -c $MODULE_FILES
ALL_OBJECTS+=$(basename --multiple $MODULE_FILES | sed 's/\.f90/\.o/g')

echo "Compiling other sources..."
$FC -c $SOURCE_FILES
ALL_OBJECTS+=" "$(basename --multiple $SOURCE_FILES | sed 's/\.f90/\.o/g')

echo "Linking..."
$FC -o aomc $ALL_OBJECTS

echo "Cleaning up intermediate files..."
rm *.o *.mod

echo "Build complete! Executable 'aomc' is ready in the root directory."
