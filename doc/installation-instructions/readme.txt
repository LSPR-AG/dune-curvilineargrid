---------------------------------------------------
-      Dune-CurvilinearGrid  v1.0                 -
---------------------------------------------------
-      Copyright                                  -
-      Aleksejs Fomins (aleksejs.fomins@lspr.ch)  -
-      Benedikt Oswald (benedikt.oswald@lspr.ch)  -
-                                                 -
-      LSPR AG                                    -
-      Grubenstrasse 9, 8045 Zurich, Switzerland  -
-      Phone: +41 43 366 90 74                    -
-                                                 -        
-      All rights reserved                        -
-                                                 -
---------------------------------------------------




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intrallation Instructions: %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

0. Introduction
1. Prerequisites
2. Installation procedure

Appendix
a1. Example .opts file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Introduction            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dune-CurvilinearGrid is a self-consistent grid manager, developed within Dune (www.dune-project.org) environment.
It is currently designed to model parallel 3D curvilinear tetrahedral geometries. Please refer to the user manual
for complete information on Dune-CurvilinearGrid. Dune-CurvilinearGeometry is a pre-requisite of
Dune-CurvilinearGrid, implementing geometric properties of 1D, 2D and 3D curvilinear simplices (e.g.
interpolation, mapping, integration). Currently, the modules can be checked out from github repository:

git clone https://github.com/LSPR-AG/dune-curvilineargeometry.git
git clone https://github.com/LSPR-AG/dune-curvilineargrid.git

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Prerequisites           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dune-CurvilinearGrid has been tested on Linux (Ubuntu 14.04, Debian 7) and OSX(10.8, 10.9, 10.10, 10.11) using
gcc 4.8, 4.9, 5.0, 5.1 and 5.2.

Current dune-curvilineargrid pre-requisites include:
1) C++ compiler supporting std=c++11 code standard
2) OpenMPI compiler
3) CMake ver >= 3.0.0 (automake/autoconf not supported)
4) Metis 5.1.0 as external library  (http://glaros.dtc.umn.edu/gkhome/metis/metis/download)
5) Parmetis 4.0.3 as external library  (http://glaros.dtc.umn.edu/gkhome/metis/parmetis/download)
Additional requirement currently is to place metis.h into include/ directory of parmetis, next to parmetis.h.
This is a quick-fix of a cmake routine that should be fixed in one of the upcoming releases.

dune-curvilienargrid also depends on Dune core modules (https://dune-project.org/downloadgit.html)
of version 2.4 (earlier versions not supported). Necessary modules include

     dune-common
     dune-geometry
     dune-grid
     dune-curvilineargeometry
     dune-curvilineargrid


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Installation procedure  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

1) Install the above external libraries
2) Ensure that paths and library paths to your external libraries are provided in .bashrc via $PATH and
    $LD_LIBRARY_PATH
2) Clone the above internal libraries into a new folder
3) Create mynicename.opts in that folder (example given below)
4) Run ./dune-common/bin/dunecontrol --use-cmake --opts=mynicename.opts all








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a1. Example .opts file     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


USE_CMAKE=yes

CONFIGURE_FLAGS="--enable-parallel --enable-experimental-grid-extensions CC=mpicc CXX=mpicxx" 

CMAKE_PREFIX_PATH="\
/opt/extlib/parmetis/4.0.3/openmpi/1.10.0/gcc/5.2.0;\
/opt/extlib/metis/5.1.0/gcc/5.2.0;\
"

## Debug opts
## GXX_WARNING_OPTS="-Wall -pedantic" 
## GXX_OPTS="-O0 -g3"

## Release opts
GXX_WARNING_OPTS=""
GXX_OPTS="-O2"

CMAKE_FLAGS=" \
-DCMAKE_CXX_FLAGS=\"$GXX_WARNING_OPTS $GXX_OPTS \" \
-DHADES_FLAGS=\" -DHAVE_DEBUG -DHAVE_PARMETIS -DHAVE_CURVGRID\" \
-DMETIS_ROOT=\"/opt/extlib/metis/5.1.0/gcc/5.2.0\" \
-DPARMETIS_ROOT=\"/opt/extlib/parmetis/4.0.3/openmpi/1.10.0/gcc/5.2.0\" \
-DCMAKE_PREFIX_PATH=\"$CMAKE_PREFIX_PATH\" \
-DDUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS:BOOL=TRUE \
-DCMAKE_SHARED_LINKER_FLAGS=\"-lzlib\" \
" 







