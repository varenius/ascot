# Install all required libraries (versions noted are tested, others may work):
- gfortran (7.5.0)
- lapack (sudo apt-get install liblapack-dev)
- OpenBLAS (sudo apt-get install libopenblas-dev)
- qt5 (5.9.5)
- hdf5 (1.10.0-patch1+docs-4)
- curl (7.58.0)
- libconfig++ (1.5-0.4; sudo apt-get install libconfig++-dev)
- boost (sudo apt-get install libboost-dev libboost-regex-dev)
- matio (sudo apt-get install libmatio-dev)
- netcdf (The old version which has c++ headers is needed: sudo apt-get install libnetcdf-cxx-legacy-dev)
- cgal (4.11-2build1; sudo apt-get install libcgal-dev)
- gurobi (8.0 or 9.1.1 seems to work both)
- libfftw3-dev (3.3.7-1)
- tclap (sudo apt-get install libtclap-dev)
- doxygen (sudo apt-get install doxygen)

# Set variables 
- $OB_DIR to where OpenBLAS is installed, if not in standard location
- $ASCOT_DIR to folder where you want to install ASCOT (e.g. /opt/ascot)
- $CPATH and $LD_LIBRARY_PATH contains paths to all relevant header files and libraries
- Edit Makefiles in subdirectories of src, and make_ascot, to correct paths etc.
- Can find qt-header locations by e.g. "qmake -query QT_INSTALL_HEADERS"
  and then set as "export CPATH=/usr/include/x86_64-linux-gnu/qt5"
- May need to "cd src/qt_plot && rm .qmake.stash && qmake" to refresh qmake stash

# Run "make_ascot all", alternatively "make_ascot all n" to run n paralell jobs - might work, if so it's faster.

# Run "bin/get_external_data" to download VMF1 etc.

Troubleshooting:
- Errors like "undefined reference to `GRBModel::read" for e.g. ivg_indep can be solved by following instructions at 
  https://support.gurobi.com/hc/en-us/articles/360039093112-How-do-I-resolve-undefined-reference-errors-while-linking-Gurobi-in-C-
  i.e. running "cd /opt/gurobi911/linux64/src/build && make && cp libgurobi_c++.a ../../lib/"
- If Segmentation Fault occurs when running ascot, please double check input config file very carefully.
- Note that "ivssrc" reference frame appears always required to be defined e.g. as 
  ("ivssrc"     , "IVSSRC", "/opt/ascot/apriori_files/IVS_SrcNamesTable.txt"),
  in your config.
