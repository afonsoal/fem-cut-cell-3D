 # This file is part of the fem-cut-cell-3D software, built using the deal.ii
 # library. You are free to use it under the GNU Lesser General Public License
 # as described in the LICENSE File.
 # Copyright (c) <2015> <Afonso Alborghetti Londero>

##
#  CMake script for the 3D Poisson problem tutorial program:
##

# Set the name of the project and target:
SET(TARGET "poisson_problem")

# Declare all source files the target consists of:
SET(TARGET_SRC
  ${TARGET}.cc

../../include/NewMesh_3D.cpp
../../include/NewCell_3D.cpp
../../include/NewFace_3D.cpp
../../include/Write_VTK.cpp
../../include/CutCell_Integration_3D.cpp

  # You can specify additional files here!
  )

# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()
