# ASCOT
ASCOT - Analysis, Scheduling and Combination Toolbox, for VLBI. Originally developed as ivg::ASCOT by the VLBI group of the Institute of Geodesy and Geoinformation of the University of Bonn (IGG). As of year 2020, ASCOT development and maintenance has been moved to the group for Space Geodesy and Geodynamics, Onsala Space Observatory, Sweden.

## Publications
- Artz et al. (2016): https://ui.adsabs.harvard.edu/abs/2016ivs..conf..217A/abstract

## Historical background
The VLBI group of the Institute of Geodesy and Geoinformation of the University of Bonn (IGG) started implementing a new analysis toolbox for VLBI observations. The main reason is the need for a flexible environment, which allows for straightforward implementations of new scientific and software-related ideas for VLBI data analysis. Furthermore, we wanted to accumulate the developments, which have been performed in Bonn in recent years, under a unified software package. The software is implemented in C++ and should finally be able to perform schedules of VLBI sessions, simulation of VLBI observations as well as geodetic data analysis and intra-technique combination. Thus, it is named: IGG VLBI Group – Analysis, Scheduling and Combination Toolbox (ivg::ASCOT). Currently, we are able to perform single-session data analysis, at a stage where the ambiguities have been resolved. Furthermore, global solutions to derive celestial and terrestrial reference frames can be performed on the normal equation level. Intra-technique combinations of several solutions complete the initial functionality of the software package. Scheduling can be performed for INTensives and small regional networks, using different approaches like impact factors or sky-coverage. 

## Install instructions
ASCOT is somewhat tricky to compile. Hopefully this will be improved in the near future, but meanwhile here's some snippets of useful information.

### Install all required libraries (versions noted are tested, others may work):
- gfortran (7.5.0)
- liblapack-dev (3.9.0-1build1)
- libopenblas-dev (0.3.8+ds-1ubuntu0.20.04.1)
- qtbase5-dev (5.9.5)
- libhdf5-dev (1.10.0-patch1+docs-4)
- curl (7.58.0)
- libconfig++-dev (1.5-0.4build1)
- libboost-dev (1.71.0-6ubuntu6)
- libboost-regex-dev (1.71.0-6ubuntu6)
- libmatio-dev (1.5.17-3)
- libnetcdf-cxx-legacy-dev (4.2-11build2; This old version which has c++ headers is needed)
- libcgal-dev (5.0.2-3)
- gurobi (8.0 or 9.1.1 seems to work both)
- libfftw3-dev (3.3.7-1)
- libtclap-dev (1.2.2-1)
- doxygen (1.8.17-0ubuntu2)
- libgsl-dev (2.5+dfsg-6build1)

### Set variables 
- $OB_DIR to where OpenBLAS is installed, if not in standard location
- $ASCOT_DIR to folder where you want to install ASCOT (e.g. /opt/ascot)
- $CPATH and $LD_LIBRARY_PATH contains paths to all relevant header files and libraries
- Edit Makefiles in subdirectories of src, and make_ascot, to correct paths etc.
- Can find qt-header locations by e.g. "qmake -query QT_INSTALL_HEADERS"
  and then set as "export CPATH=/usr/include/x86_64-linux-gnu/qt5"
- Need to "cd src/qt_plot && rm .qmake.stash && qmake" to refresh qmake stash
- Don't forget to add gurobi to your path, e.g.
  export GUROBI_HOME="/opt/gurobi911/linux64"
  export PATH="${PATH}:${GUROBI_HOME}/bin"
  export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
  export GRB_LICENSE_FILE="/opt/gurobi911/gurobi.lic"

### Run "make_ascot all", alternatively "make_ascot all n" to run n paralell jobs - might work, if so it's faster.

### Run "bin/get_external_data" to download VMF1 etc.

## Troubleshooting
- Errors like "undefined reference to `GRBModel::read" for e.g. ivg_indep can be solved by following instructions at 
  https://support.gurobi.com/hc/en-us/articles/360039093112-How-do-I-resolve-undefined-reference-errors-while-linking-Gurobi-in-C-
  i.e. running "cd /opt/gurobi911/linux64/src/build && make && cp libgurobi_c++.a ../../lib/"
- If Segmentation Fault occurs when running ascot, please double check input config file very carefully.
- Note that "ivssrc" reference frame appears always required to be defined e.g. as 
  ("ivssrc"     , "IVSSRC", "/opt/ascot/apriori_files/IVS_SrcNamesTable.txt"), in your config.
- If ivg compilation fails with QtCore/QString problems, remember to "cd src/qt_plot && rm .qmake.stash && qmake" to refresh qmake stash.
