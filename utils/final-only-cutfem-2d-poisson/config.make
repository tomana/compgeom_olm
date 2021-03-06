################################################################################
# CONFIGURE PROJECT MAKEFILE (optional)
#   This file is where we make project specific configurations.
################################################################################

################################################################################
# OF ROOT
#   The location of your root openFrameworks installation
#       (default) OF_ROOT = ../../../../openFrameworks/../openFrameworks/../openFrameworks/../openFrameworks/../openFrameworks/../openFrameworks/../openFrameworks/../openFrameworks/../openFrameworks/../openFrameworks/../openFrameworks/../../openFrameworks/../../openFrameworks/../../openFrameworks/../../openFrameworks/../../openFrameworks/../../openFrameworks
################################################################################
OF_ROOT = /home/tom/Desktop/openFrameworks

################################################################################
# PROJECT ROOT
#   The location of the project - a starting place for searching for files
#       (default) PROJECT_ROOT = . (this directory)
#
################################################################################
# PROJECT_ROOT = .

################################################################################
# PROJECT SPECIFIC CHECKS
#   This is a project defined section to create internal makefile flags to
#   conditionally enable or disable the addition of various features within
#   this makefile.  For instance, if you want to make changes based on whether
#   GTK is installed, one might test that here and create a variable to check.
################################################################################
# None

################################################################################
# PROJECT EXTERNAL SOURCE PATHS
#   These are fully qualified paths that are not within the PROJECT_ROOT folder.
#   Like source folders in the PROJECT_ROOT, these paths are subject to
#   exlclusion via the PROJECT_EXLCUSIONS list.
#
#     (default) PROJECT_EXTERNAL_SOURCE_PATHS = (blank)
#
#   Note: Leave a leading space when adding list items with the += operator
################################################################################
#PROJECT_EXTERNAL_SOURCE_PATHS =

################################################################################
# PROJECT EXCLUSIONS
#   These makefiles assume that all folders in your current project directory
#   and any listed in the PROJECT_EXTERNAL_SOURCH_PATHS are are valid locations
#   to look for source code. The any folders or files that match any of the
#   items in the PROJECT_EXCLUSIONS list below will be ignored.
#
#   Each item in the PROJECT_EXCLUSIONS list will be treated as a complete
#   string unless teh user adds a wildcard (%) operator to match subdirectories.
#   GNU make only allows one wildcard for matching.  The second wildcard (%) is
#   treated literally.
#
#      (default) PROJECT_EXCLUSIONS = (blank)
#
#		Will automatically exclude the following:
#
#			$(PROJECT_ROOT)/bin%
#			$(PROJECT_ROOT)/obj%
#			$(PROJECT_ROOT)/%.xcodeproj
#
#   Note: Leave a leading space when adding list items with the += operator
################################################################################
# PROJECT_EXCLUSIONS =

################################################################################
# PROJECT LINKER FLAGS
#	These flags will be sent to the linker when compiling the executable.
#
#		(default) PROJECT_LDFLAGS = -Wl,-rpath=./libs
#
#   Note: Leave a leading space when adding list items with the += operator
################################################################################

# Currently, shared libraries that are needed are copied to the
# $(PROJECT_ROOT)/bin/libs directory.  The following LDFLAGS tell the linker to
# add a runtime path to search for those shared libraries, since they aren't
# incorporated directly into the final executable application binary.
# TODO: should this be a default setting?
PROJECT_LDFLAGS= -lgl2ps -std=c++11 -pg -g   -std=c++11  -Wno-deprecated -frounding-math  -fopenmp -O2 -g -DNDEBUG -rdynamic /home/tom/Desktop/FEniCS-cutfem/cutfem/local.naerland.overlappingmeshes/lib/libcutfem.so /home/tom/Desktop/FEniCS-cutfem/dolfin-cutfem/local.massing.topic-cutfem/lib/libdolfin.so -lxml2 -lboost_filesystem-mt -lboost_program_options-mt -lboost_system-mt -lboost_thread-mt -lpthread -lboost_iostreams-mt -lhdf5 -lpthread -lz -lrt -lm /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libml.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libgaleri-xpetra.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libgaleri.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libisorropia.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libxpetra-sup.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libxpetra-ext.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libxpetra.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libzoltan.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libifpack.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libaztecoo.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libamesos.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libepetraext.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libtriutils.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libbelostpetra.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libbelosepetra.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libbelos.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libtpetraext.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libtpetrainout.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libtpetra.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libepetra.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libkokkosdisttsqr.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libkokkosnodetsqr.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libkokkoslinalg.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libkokkosnodeapi.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libkokkos.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libtpi.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libteuchosremainder.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libteuchosnumerics.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libteuchoscomm.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libteuchosparameterlist.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libteuchoscore.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libslepc.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libpetsc.so -lumfpack -lamd -lcblas -lf77blas -latlas -lcholmod -lamd -lcamd -lcolamd -lccolamd -lrt /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libparmetis.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libmetis.so -llapack -lcblas -lf77blas -latlas -lcblas -lf77blas -latlas -lgfortran -lgfortran -lcholmod -lamd -lcamd -lcolamd -lccolamd -lrt /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libparmetis.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libmetis.so -llapack -lcblas -lf77blas -latlas -lcblas -lf77blas -latlas -lgfortran /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libptscotch.a /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libscotch.a /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libptscotcherr.a /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libparmetis.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libmetis.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libCGAL.so -lboost_thread-mt -lpthread -lboost_system-mt -lgmp -lmpfr -lz -lcppunit /usr/lib/openmpi/lib/libmpi_cxx.so /usr/lib/openmpi/lib/libmpi.so /usr/lib/openmpi/lib/libopen-rte.so /usr/lib/openmpi/lib/libopen-pal.so -ldl -lnsl -lutil -lm -ldl -lQtGui -lQtCore -lboost_system-mt -lboost_thread-mt -lpthread -lboost_iostreams-mt -lhdf5 -lpthread -lz -lrt /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libml.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libgaleri-xpetra.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libgaleri.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libisorropia.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libxpetra-sup.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libxpetra-ext.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libxpetra.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libzoltan.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libifpack.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libaztecoo.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libamesos.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libepetraext.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libtriutils.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libbelostpetra.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libbelosepetra.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libbelos.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libtpetraext.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libtpetrainout.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libtpetra.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libepetra.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libkokkosdisttsqr.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libkokkosnodetsqr.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libkokkoslinalg.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libkokkosnodeapi.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libkokkos.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libtpi.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libteuchosremainder.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libteuchosnumerics.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libteuchoscomm.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libteuchosparameterlist.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libteuchoscore.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libslepc.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libpetsc.so -lumfpack -lamd -lcblas -lf77blas -latlas -lcholmod -lcamd -lcolamd -lccolamd /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libparmetis.so /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libmetis.so -llapack -lgfortran /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libptscotch.a /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libscotch.a /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libptscotcherr.a /home/tom/Desktop/FEniCS-dev/FEniCS/lib/libCGAL.so -lgmp -lmpfr -lcppunit /usr/lib/openmpi/lib/libmpi_cxx.so /usr/lib/openmpi/lib/libmpi.so /usr/lib/openmpi/lib/libopen-rte.so /usr/lib/openmpi/lib/libopen-pal.so -lnsl -lutil -lQtGui -lQtCore -Wl,-rpath,/home/tom/Desktop/FEniCS-cutfem/cutfem/local.naerland.overlappingmeshes/lib:/home/tom/Desktop/FEniCS-cutfem/dolfin-cutfem/local.massing.topic-cutfem/lib:/home/tom/Desktop/FEniCS-dev/FEniCS/lib:/usr/lib/openmpi/lib
################################################################################
# PROJECT DEFINES
#   Create a space-delimited list of DEFINES. The list will be converted into
#   CFLAGS with the "-D" flag later in the makefile.
#
#		(default) PROJECT_DEFINES = (blank)
#
#   Note: Leave a leading space when adding list items with the += operator
################################################################################
PROJECT_DEFINES =
-DBOOST_UBLAS_NDEBUG -DCUTFEM_VERSION=\"0.0.1+\" -DDEBUG -DENABLE_PETSC_SNES -DHAS_CGAL -DHAS_CHOLMOD -DHAS_CPPUNIT -DHAS_HDF5 -DHAS_MPI -DHAS_OPENMP -DHAS_PARMETIS -DHAS_PETSC -DHAS_QT4 -DHAS_SCOTCH -DHAS_SLEPC -DHAS_TRILINOS -DHAS_UMFPACK -DHAS_VTK -DHAS_ZLIB -D_BSD_SOURCE -D_FORTIFY_SOURCE=2 -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE

################################################################################
# PROJECT CFLAGS
#   This is a list of fully qualified CFLAGS required when compiling for this
#   project.  These CFLAGS will be used IN ADDITION TO the PLATFORM_CFLAGS
#   defined in your platform specific core configuration files. These flags are
#   presented to the compiler BEFORE the PROJECT_OPTIMIZATION_CFLAGS below.
#
#		(default) PROJECT_CFLAGS = (blank)
#
#   Note: Before adding PROJECT_CFLAGS, note that the PLATFORM_CFLAGS defined in
#   your platform specific configuration file will be applied by default and
#   further flags here may not be needed.
#
#   Note: Leave a leading space when adding list items with the += operator
################################################################################
PROJECT_CFLAGS = -std=c++11 -pg -g   -std=c++11  -Wno-deprecated   -frounding-math  -fopenmp -O2 -g -DNDEBUG -I/home/tom/Desktop/FEniCS-cutfem/cutfem/local.naerland.overlappingmeshes/include -I/home/tom/Desktop/FEniCS-cutfem/dolfin-cutfem/local.massing.topic-cutfem/include -I/usr/include/libxml2 -I/home/tom/Desktop/FEniCS-cutfem/ffc-cutfem/build/include -I/home/tom/Desktop/FEniCS-dev/FEniCS -I/home/tom/Desktop/FEniCS-dev/FEniCS/include -I/usr/include/suitesparse -I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi -I/usr/include/eigen3 -I/home/tom/Desktop/FEniCS-dev/FEniCS/lib/cmake/Trilinos/../../../../openFrameworks/../openFrameworks/../openFrameworks/include/trilinos -I/home/tom/Desktop/FEniCS-dev/FEniCS/include/trilinos -I/usr/include/qt4 -I/usr/local/include/vtk-6.1     -DDOLFIN_DEPRECATION_ERROR

#PROJECT_CFLAGS = -DCUTFEM_VERSION=\"0.0.1+\" -DHAS_CGAL -DHAS_CHOLMOD -DHAS_CPPUNIT -DHAS_MPI -DHAS_OPENMP -DHAS_PARMETIS -DHAS_PETSC -DHAS_SCOTCH -DHAS_SLEPC -DHAS_TRILINOS
#-std=c++11 -g -Wno-deprecated   -frounding-math  -fopenmp -I/home/tom/Desktop/FEniCS-cutfem/cutfem/local.naerland.olm-geometry/include
#-I/home/tom/Desktop/FEniCS-cutfem/dolfin-cutfem/local.massing.topic-cutfem/include -I/usr/include/libxml2 -I/home/tom/Desktop/FEniCS-cutfem/ufc-cutfem/build/include
#-I/home/tom/Desktop/FEniCS-dev/FEniCS -I/home/tom/Desktop/FEniCS-dev/FEniCS/include -I/usr/include/suitesparse -I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi
#-I/usr/include/eigen3 -I/home/tom/Desktop/FEniCS-dev/FEniCS/lib/cmake/Trilinos/../../../../openFrameworks/../openFrameworks/../openFrameworks/../openFrameworks/../openFrameworks/include/trilinos -I/home/tom/Desktop/FEniCS-dev/FEniCS/include/trilinos
#-DDOLFIN_DEPRECATION_ERROR
################################################################################
# PROJECT OPTIMIZATION CFLAGS
#   These are lists of CFLAGS that are target-specific.  While any flags could
#   be conditionally added, they are usually limited to optimization flags.
#   These flags are added BEFORE the PROJECT_CFLAGS.
#
#   PROJECT_OPTIMIZATION_CFLAGS_RELEASE flags are only applied to RELEASE targets.
#
#		(default) PROJECT_OPTIMIZATION_CFLAGS_RELEASE = (blank)
#
#   PROJECT_OPTIMIZATION_CFLAGS_DEBUG flags are only applied to DEBUG targets.
#
#		(default) PROJECT_OPTIMIZATION_CFLAGS_DEBUG = (blank)
#
#   Note: Before adding PROJECT_OPTIMIZATION_CFLAGS, please note that the
#   PLATFORM_OPTIMIZATION_CFLAGS defined in your platform specific configuration
#   file will be applied by default and further optimization flags here may not
#   be needed.
#
#   Note: Leave a leading space when adding list items with the += operator
################################################################################
# PROJECT_OPTIMIZATION_CFLAGS_RELEASE =
# PROJECT_OPTIMIZATION_CFLAGS_DEBUG =

################################################################################
# PROJECT COMPILERS
#   Custom compilers can be set for CC and CXX
#		(default) PROJECT_CXX = (blank)
#		(default) PROJECT_CC = (blank)
#   Note: Leave a leading space when adding list items with the += operator
################################################################################
# PROJECT_CXX =
# PROJECT_CC =
