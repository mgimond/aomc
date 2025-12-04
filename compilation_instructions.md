# Compiling the AOMC Fortran Code

This document provides instructions on how to compile the Fortran source code for the AOMC project using the `gfortran` compiler.

## 1. Prerequisites

You must have the GNU Fortran compiler (`gfortran`) installed on your system. You can check if it's installed by running:

```sh
gfortran --version
```

## 2. Understanding Dependencies

The source code uses Fortran `MODULE`s to share procedures and data. Files that define modules must be compiled before the files that `USE` them.

- **Module Definition Files:** `global.f90`, `modules.f90`

These files will be compiled first to generate the necessary `.mod` files, which describe the module interfaces to the compiler.

## 3. Compilation Steps

The compilation process involves two main stages:
1.  Compiling individual source files (`.f90`) into object files (`.o`).
2.  Linking the object files together to create the final executable.

### Step 3.1: Compile Module Files

First, compile the files that define the modules. This will create `global.o`, `modules.o`, and the corresponding `.mod` files.

```sh
gfortran -c global.f90 modules.f90
```

### Step 3.2: Compile Remaining Source Files

Next, compile the rest of the source files. The compiler will use the `.mod` files created in the previous step to resolve dependencies.

```sh
gfortran -c geom2.f90 interface.f90 light_internal.f90 logbin.f90 mc.f90 vsf.f90 water.f90
```

### Step 3.3: Link Object Files

Finally, link all the generated object files (`.o`) together to create the executable named `aomc`.

```sh
gfortran -o aomc *.o
```

### Step 3.4: Clean Up (Optional)

You can remove the intermediate object (`.o`) and module (`.mod`) files after the executable is created.

```sh
rm *.o *.mod
```

## 4. Automation Script

To simplify the process, you can use the following shell script. Save it as `build.sh` and make it executable (`chmod +x build.sh`).

```sh
#!/bin/bash
# A script to compile the AOMC Fortran project.

# Exit immediately if a command exits with a non-zero status.
set -e

# Compiler
FC=gfortran

# List of source files
MODULE_FILES="global.f90 modules.f90"
SOURCE_FILES="geom2.f90 interface.f90 light_internal.f90 logbin.f90 mc.f90 vsf.f90 water.f90"
ALL_OBJECTS=""

echo "Compiling modules..."
$FC -c $MODULE_FILES
ALL_OBJECTS+=$(echo $MODULE_FILES | sed 's/\.f90/\.o/g')

echo "Compiling other sources..."
$FC -c $SOURCE_FILES
ALL_OBJECTS+=" "$(echo $SOURCE_FILES | sed 's/\.f90/\.o/g')

echo "Linking..."
$FC -o aomc $ALL_OBJECTS

echo "Cleaning up intermediate files..."
rm *.o *.mod

echo "Build complete! Executable 'aomc' is ready."
```
