project(poc-solvers)
cmake_minimum_required(VERSION 2.8)

####################################################################################################
# NOTE: 
#
#  INSTALL.sh setup and installation:
#
#  Expert  and installation:
#
#    set adequately LD_LIBRARY_PATH, CMAKE_LIBRARY_PATH and CMAKE_INCLUDE_PATH environment variables
#    so locally installed libraries and includes can be found
#
#    MUMPSLIB=/home/softs/MUMPS_5.1.0/lib
#    MUMPSINC=/home/softs/MUMPS_5.1.0/include
#
#    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MUMPSLIB
#
#    export CMAKE_LIBRARY_PATH=$LD_LIBRARY_PATH
#    export CMAKE_INCLUDE_PATH=$CMAKE_INCLUDE_PATH:$MUMPSINC
#
#    warning : LD_LIBRARY_PATH and CMAKE_LIBRARY_PATH should contain the same paths.
#
#
#
#  mkdir cmade; cd cmade; cmake -DCMAKE_BUILD_TYPE=Debug ../src && make install
#
#  -DLOCAL_INSTALL=<PATH> where PATH is root path to companion softs installation directory
#
#  Examples:
#
#  cmake ../src -DLOCAL_INSTALL=/home/softs/local/ -DOLD_ZDHIPS=1 -DMUMPS=NO -DPASTIX=NO -DHIPS=NO -DUMFPACK=NO
#
#
####################################################################################################

set(STATIC NO CACHE BOOL "Set to YES to have a static library. Disables PASTIX.")
set(HIPS YES CACHE BOOL "Set to YES to auto-detect HIPS. If found, it will enable MPI.")
set(ARPACK YES CACHE BOOL "Set to YES to auto-detect ARPACK.")
set(PARPACK YES CACHE BOOL "Set to YES to auto-detect PARPACK. If found, it will enable MPI.")
set(PASTIX YES CACHE BOOL "Set to YES to enable PASTIX. If it can not find it or if PASTIX is disabled, it will enable UMFPACK instead. PASTIX may enable MPI.")
set(MUMPS YES CACHE BOOL "Set to YES to enable MUMPS. It will enable MPI and SCALAPACK.")
set(MPI NO CACHE BOOL "Set to YES to enable MPI. Enabled by MUMPS, PARPACK or HIPS anyway.")
set(UMFPACK YES CACHE BOOL "Set to YES to enable UMFPACK. Needed by TUGO. It will enable LAPACK.")
set(SCALAPACK NO CACHE BOOL "Set to YES to enable SCALAPACK. Only needed and enabled by MUMPS anyway. It will enable LAPACK and MPI.")
set(LAPACK YES CACHE BOOL "Set to YES to enable LAPACK. Enabled by UMFPACK and SCALAPACK anyway. It will enable BLAS.")
set(BLAS NO CACHE BOOL "Set to YES to enable BLAS. Enabled by LAPACK anyway.")
set(PETSC YES CACHE BOOL "Set to YES to enable PETSC")

if(LOCAL_INSTALL)
else()
  execute_process(COMMAND readlink -f ../.. OUTPUT_VARIABLE LOCAL_INSTALL OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()
message("------------> LOCAL_INSTALL : ${LOCAL_INSTALL}")

################################################################################
# CXSparse

# make -C ../CXSparse/Lib distclean # if you have suitesparse-devel or libsuitesparse-dev installed
find_library(CSPARSELIB NAMES cxsparse HINTS ../CXSparse/Lib)

if(CSPARSELIB)
  message("+++ CSPARSELIB : ${CSPARSELIB}")
else()
  message("!!! compiling CSPARSE library")
  execute_process(COMMAND make -C ../CXSparse/Lib RESULT_VARIABLE CSPARSELIBMADE)
  message("+++ ${CSPARSELIBMADE}")
  
  if(CSPARSELIBMADE EQUAL 0)
    message(FATAL_ERROR "Run the same cmake command again and it should work.")
  else()
    message(FATAL_ERROR "specify option `-DCSPARSELIB=../path/to/libcxsparse.a'")
  endif()

endif()

################################################################################
# PaStiX

if(PASTIX AND NOT STATIC)
  find_program(PASTIXCONF pastix-conf HINTS ${LOCAL_INSTALL}/bin)
  message("+++ PASTIXCONF : ${PASTIXCONF}")
  
  if(PASTIXCONF)
    execute_process(COMMAND ${PASTIXCONF} --incs OUTPUT_VARIABLE PASTIX_CPP)
    add_definitions(-DPASTIX ${PASTIX_CPP})
    message("+++ PASTIXINCS : ${PASTIX_CPP}")
    # depended upon libraries
    execute_process(COMMAND ${PASTIXCONF} --libs OUTPUT_VARIABLE PASTIXLIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
    find_package(OpenMP) # /usr/share/doc/packages/cmake/cmake.html#module:FindOpenMP
    message(WARNING "!!! prepending `-lblas ${OpenMP_C_FLAGS}' to PASTIXLIBS !!!")
    set(PASTIXLIBS -lblas ${OpenMP_C_FLAGS} ${PASTIXLIBS})
    message("+++ PASTIXLIBS : ${PASTIXLIBS}")
    message(WARNING "!!! prepending `-lmatrix_driver' to TESTLIBS !!!")
    set(TESTLIBS -lmatrix_driver ${TESTLIBS})
    message("+++ TESTLIBS : ${TESTLIBS}")
    # check if MPI is needed
    execute_process(COMMAND pastix-conf --cc OUTPUT_VARIABLE PASTIXCC)
    
    if(${PASTIXCC} MATCHES mpi)
      message("+++ Your PASTIX version needs MPI")
      set(MPI YES)
      message(WARNING "!!! appending `-open-pal' to PASTIXLIBS !!!")
      set(PASTIXLIBS ${PASTIXLIBS} -lopen-pal)
    endif()
    
    message(STATUS "*** PASTIX library found")
  else()
    message(WARNING "*** PASTIX library not found : disabling.")
    set(PASTIX NO)
  endif()
  
endif()

#message(STATUS "MPI= ${MPI}")

if(NOT PASTIX AND NOT UMFPACK)
  message(WARNING "*** PASTIX disabled : enabling UMFPACK instead")
  set(UMFPACK YES)
endif()

################################################################################
# HIPS

# if(HIPS)
#   find_library(HIPSLIB NAMES hips HINTS ${LOCAL_INSTALL}/hips/LIB)
#   find_library(HIPS_io_LIB NAMES io HINTS ${LOCAL_INSTALL}/hips/LIB)
#   find_library(HIPS_spkit_LIB NAMES spkit HINTS ${LOCAL_INSTALL}/hips/LIB)
#   
#   message(STATUS "found : ${HIPSLIB} ${HIPS_io_LIB} ${HIPS_spkit_LIB}")
#   if(HIPSLIB AND HIPS_io_LIB AND HIPS_spkit_LIB)
#     add_definitions(-DHIPS)
#     # add_definitions(-D__DEBUG_HIPS)
#     set(HIPSLIBS ${HIPSLIB} ${HIPS_io_LIB} ${HIPS_spkit_LIB})
#     find_path(HIPSINCDIRS NAMES hips.h PATHS ${LOCAL_INSTALL}/hips/LIB)
#     include_directories(${HIPSINCDIRS})
#     message("+++ HIPS needs MPI")
#     set(MPI YES)
#   else()
#     message(WARNING "*** HIPS library not found")
#   endif()
#   
# endif()

################################################################################
# ZDHIPS

if(HIPS)
  
  if(OLD_ZDHIPS)
    message(STATUS "testing : ${LOCAL_INSTALL}/zdhips/dhips/DLIB ${LOCAL_INSTALL}/zdhips/zhips/ZLIB")
    
    find_library(DHIPSLIB NAMES dhips HINTS ${LOCAL_INSTALL}/zdhips/dhips/DLIB)
    find_library(DHIPS_io_LIB NAMES dio HINTS ${LOCAL_INSTALL}/zdhips/dhips/DLIB)
    find_library(DHIPS_spkit_LIB NAMES dspkit HINTS ${LOCAL_INSTALL}/zdhips/dhips/DLIB)
    
    find_library(ZHIPSLIB NAMES zhips HINTS ${LOCAL_INSTALL}/zdhips/zhips/ZLIB)
    find_library(ZHIPS_io_LIB NAMES zio HINTS ${LOCAL_INSTALL}/zdhips/zhips/ZLIB)
    find_library(ZHIPS_spkit_LIB NAMES zspkit HINTS ${LOCAL_INSTALL}/zdhips/zhips/ZLIB)
  else()
    message(STATUS "testing : ${LOCAL_INSTALL}/zdhips/dhips/LIB ${LOCAL_INSTALL}/zdhips/zhips/LIB")
    
    find_library(DHIPSLIB NAMES hips HINTS ${LOCAL_INSTALL}/zdhips/dhips/LIB)
    find_library(DHIPS_io_LIB NAMES io HINTS ${LOCAL_INSTALL}/zdhips/dhips/LIB)
    find_library(DHIPS_spkit_LIB NAMES spkit HINTS ${LOCAL_INSTALL}/zdhips/dhips/LIB)
    
    find_library(ZHIPSLIB NAMES hips HINTS ${LOCAL_INSTALL}/zdhips/zhips/LIB)
    find_library(ZHIPS_io_LIB NAMES io HINTS ${LOCAL_INSTALL}/zdhips/zhips/LIB)
    find_library(ZHIPS_spkit_LIB NAMES spkit HINTS ${LOCAL_INSTALL}/zdhips/zhips/LIB)
  endif()
  
  message(STATUS "found : ${DHIPSLIB} ${DHIPS_io_LIB} ${DHIPS_spkit_LIB} ${ZHIPSLIB} ${ZHIPS_io_LIB} ${ZHIPS_spkit_LIB}")
  if(DHIPSLIB AND DHIPS_io_LIB AND DHIPS_spkit_LIB AND ZHIPSLIB AND ZHIPS_io_LIB AND ZHIPS_spkit_LIB)
    add_definitions(-DHIPS)
    # add_definitions(-D__DEBUG_HIPS)
    set(HIPSLIBS ${DHIPSLIB} ${DHIPS_io_LIB} ${DHIPS_spkit_LIB} ${ZHIPSLIB} ${ZHIPS_io_LIB} ${ZHIPS_spkit_LIB})
    
    if(OLD_ZDHIPS)
      find_path(DHIPSINCDIRS NAMES dhips.h PATHS ${LOCAL_INSTALL}/zdhips/dhips/DLIB)
      find_path(ZHIPSINCDIRS NAMES zhips.h PATHS ${LOCAL_INSTALL}/zdhips/zhips/ZLIB)
    else()
      find_path(DHIPSINCDIRS NAMES dhips.h PATHS ${LOCAL_INSTALL}/zdhips/dhips/LIB)
      find_path(ZHIPSINCDIRS NAMES zhips.h PATHS ${LOCAL_INSTALL}/zdhips/zhips/LIB)
    endif()
    
    set(HIPSINCDIRS ${DHIPSINCDIRS} ${ZHIPSINCDIRS})
    include_directories(${HIPSINCDIRS})
    
    message("+++ HIPS needs MPI")
    set(MPI YES)
  else()
    message(WARNING "*** HIPS library not found")
  endif()
  
endif()

################################################################################
# ARPACK and PARPACK

if(ARPACK)
  find_library(ARPACKLIB NAMES arpack)
  
  if(ARPACKLIB)
    message(STATUS "*** ARPACK library found")
  else()
    message(WARNING "*** ARPACK library not found")
    set(ARPACKLIB "")
  endif()
  
  message("+++ ARPACKLIB : ${ARPACKLIB}")
endif()

if(PARPACK)
  find_library(PARPACKLIB NAMES parpack)
  
  if(PARPACKLIB)
    message("+++ PARPACK needs MPI")
    set(MPI YES)
    message(STATUS "*** PARPACK library found")
  else()
    message(WARNING "*** PARPACK library not found")
    set(PARPACKLIB "")
  endif()
  
  message("+++ PARPACKLIB : ${PARPACKLIB}")
endif()

################################################################################
# SCALAPACK

if(MUMPS)
  message("+++ MUMPS needs SCALAPACK")
  set(SCALAPACK YES)
endif()

if(SCALAPACK)
  find_library(SCALAPACKLIB NAMES scalapack scalapack-openmpi HINTS ${LOCAL_INSTALL}/lib)
  
  if(SCALAPACKLIB)
    message(STATUS "*** SCALAPACK library found")
    
    message("+++ SCALAPACK needs LAPACK and MPI")
    set(LAPACK YES)
    set(MPI YES)
  else()
    message(WARNING "*** SCALAPACK library not found : $(SCALAPACKLIB)")
  endif()
  
  message("+++ SCALAPACKLIB : ${SCALAPACKLIB}")
endif()


################################################################################
# MUMPS

if(MUMPS)
  # find_path(MUMPSINCDIRS NAMES dmumps_c.h PATHS /usr/include/mumps)
  # # THE FOLLOWING LINE WILL GIVE PRIORITY TO /usr/include/mumps/mpi.h
  # # WHICH WILL MESS-UP MPI !!!
  # include_directories(${MUMPSINCDIRS})
#   find_path(MUMPSINCDIRS NAMES dmumps_c.h HINTS ${LOCAL_INSTALL}/include HINTS /usr/include/mumps)
  include_directories(${LOCAL_INSTALL}/include)
#   find_library(DMUMPSLIB NAMES dmumps_seq dmumps)
#   find_library(ZMUMPSLIB NAMES zmumps_seq zmumps)
  find_library(DMUMPSLIB NAMES dmumps HINTS ${LOCAL_INSTALL}/lib HINTS /usr/lib64/mpi/gcc/openmpi/lib64/)
  find_library(ZMUMPSLIB NAMES zmumps HINTS ${LOCAL_INSTALL}/lib HINTS /usr/lib64/mpi/gcc/openmpi/lib64/)
  find_library(MUMPS_COMMONLIB NAMES mumps_common HINTS ${LOCAL_INSTALL}/lib HINTS /usr/lib64/mpi/gcc/openmpi/lib64/)
  find_library(PORDLIB NAMES pord                 HINTS ${LOCAL_INSTALL}/lib HINTS /usr/lib64/mpi/gcc/openmpi/lib64/)
  # SCOTCH libraries (yes, even esmumps) :
  find_library(ESMUMPSLIB NAMES esmumps     HINTS ${LOCAL_INSTALL}/lib)
  find_library(SCOTCHLIB NAMES scotch       HINTS ${LOCAL_INSTALL}/lib)
  find_library(PTSCOTCHLIB NAMES ptscotch   HINTS ${LOCAL_INSTALL}/lib)
  find_library(SCOTCHERRLIB NAMES scotcherr HINTS ${LOCAL_INSTALL}/lib)
  # METIS libraries :
  find_library(METISLIB NAMES metis         HINTS ${LOCAL_INSTALL}/lib)
  find_library(MPIFLIB NAMES mpi_f77 mpi_mpifh HINTS /usr/lib64/mpi/gcc/openmpi/lib64/)
  
  set(MUMPSLIBS ${DMUMPSLIB} ${ZMUMPSLIB} ${MUMPS_COMMONLIB} ${PORDLIB} ${ESMUMPSLIB} ${PTSCOTCHLIB} ${SCOTCHLIB} ${SCOTCHERRLIB} ${METISLIB} ${MPIFLIB})
#   set(MUMPSLIBS ${DMUMPSLIB} ${ZMUMPSLIB} ${MUMPS_COMMONLIB} ${PORDLIB} ${ESMUMPSLIB} ${PTSCOTCHLIB} ${SCOTCHLIB} ${SCOTCHERRLIB} ${MPIFLIB})
  if(DMUMPSLIB AND ZMUMPSLIB AND MUMPS_COMMONLIB AND PORDLIB AND ESMUMPSLIB  AND MPIFLIB)
    message(STATUS "*** MUMPS libraries found")
    add_definitions(-DMUMPS)
    
    message("+++ MUMPS needs MPI")
    set(MPI YES)
  else()
    message(WARNING "*** at least 1 MUMPS library not found : ${MUMPSLIBS}")
    set(MUMPSLIBS "")
  endif()
  
  message("+++ MUMPSLIBS : ${MUMPSLIBS}")
  add_definitions(-D__TESTLIB)
endif()

################################################################################
# LAPACK

if(UMFPACK)
  message("+++ UMFPACK may need LAPACK")
  set(LAPACK YES)
endif()

if(LAPACK)
  message("+++ LAPACK may need BLAS and GFORTRAN")
  set(GFORTRAN YES)
  set(BLAS YES)
endif()

if(GFORTRAN)
  # find_library(GFORTRANLIB NAMES gfortran) # DOES NOT WORK
  set(GFORTRANLIB -lgfortran) # WORKS
  message("+++ GFORTRANLIB : ${GFORTRANLIB}")
endif()

if(BLAS)
  find_library(BLASLIB NAMES blas)
  
  if(BLASLIB)
    message(STATUS "*** BLAS library found")
  else()
    message(WARNING "*** BLAS library not found")
  endif()
  
  message("+++ BLASLIB : ${BLASLIB}")
endif()

if(LAPACK)
  add_definitions(-DLAPACK)
  find_library(LAPACKLIB NAMES lapack)
  
  if(LAPACKLIB)
    message(STATUS "*** LAPACK library found")
  else()
    message(WARNING "*** LAPACK library not found")
  endif()
  
  message("+++ LAPACKLIB : ${LAPACKLIB}")
endif()

################################################################################
# UMFPACK

if(UMFPACK)
  add_definitions(-DUMFPACK)
  find_path(UMFPACKINCDIRS NAMES umfpack.h PATHS /usr/include/suitesparse PATHS ${LOCAL_INSTALL}/include)
  include_directories(${UMFPACKINCDIRS})
  find_library(UMFPACKLIB NAMES umfpack HINTS ${LOCAL_INSTALL}/lib)
  find_library(MLIB NAMES m)
  find_library(RTLIB NAMES rt)
  if(NOT RTLIB)
    message(WARNING "*** RT library not found: failsafing")
    set(RTLIB "")
    endif()
  set(UMFPACKLIBS ${UMFPACKLIB} ${LAPACKLIB} ${BLASLIB} ${GFORTRANLIB} ${MLIB} ${RTLIB})
  message("+++ UMFPACKLIBS : ${UMFPACKLIBS}")
  endif()

if(UMFPACK AND NOT UMFPACKLIB)
  message(WARNING "*** Can not find UMFPACK")
  endif()

if(UMFPACK AND UMFPACKLIB)
  message("+++ Auto-detecting UMFPACK dependencies")
  set(UMFPACK_TEST_FILE ${CMAKE_BINARY_DIR}/CMakeFiles/CMakeTmp/umfpacktest.c)
  set(APPENDED_SOMETHING TRUE)
  set(WORD_REGEX "[^a-zA-Z]")
  
  foreach(UMFPACK_TESTS 1 2 3)
    message(STATUS "iteration ${UMFPACK_TESTS}")
    
    message("+++ UMFPACKLIBS : ${UMFPACKLIBS}")
    execute_process(COMMAND echo "void main(){umfpack_di_solve(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);umfpack_zi_symbolic(0,0,0,0,0,0,0,0,0);umf_i_free(0);}" OUTPUT_FILE ${UMFPACK_TEST_FILE})#cleaned-up after each try_compile
    try_compile(UMFPACK_COMPILES ${CMAKE_BINARY_DIR} ${UMFPACK_TEST_FILE}
      CMAKE_FLAGS "-DLINK_LIBRARIES=${UMFPACKLIBS}"
      OUTPUT_VARIABLE UMFPACK_OUTPUT)
    message("!!! UMFPACK_COMPILES already ? - ${UMFPACK_COMPILES}")
#     message("vvvvvvvvvvvvvvvvvv--- output below ---vvvvvvvvvvvvvvvvvvvvvvvvvvv")
#     message("${UMFPACK_OUTPUT}")
#     message("^^^^^^^^^^^^^^^^^^--- output above ---^^^^^^^^^^^^^^^^^^^^^^^^^^^")
    
    if(UMFPACK_COMPILES)
      break()
    endif()
    
    set(APPENDED_SOMETHING FALSE)
    
    foreach(UMFPACK_DEPENDENCY cholmod amd camd colamd ccolamd)
#       message(STATUS "checking for ${UMFPACK_DEPENDENCY} in ${UMFPACK_TESTS}")
      
      if(UMFPACK_OUTPUT MATCHES "${WORD_REGEX}${UMFPACK_DEPENDENCY}_")
        find_library(UMFPACK_${UMFPACK_DEPENDENCY}_LIB NAMES ${UMFPACK_DEPENDENCY}  PATHS HINTS ${LOCAL_INSTALL}/lib)
        
        if(NOT UMFPACK_${UMFPACK_DEPENDENCY}_LIB)
          message(FATAL_ERROR "Can not find ${UMFPACK_DEPENDENCY} library.")
        endif()
        
        message("!!! Appending ${UMFPACK_${UMFPACK_DEPENDENCY}_LIB} for ${UMFPACK_DEPENDENCY}")
        set(UMFPACKLIBS ${UMFPACKLIBS} ${UMFPACK_${UMFPACK_DEPENDENCY}_LIB})
        set(APPENDED_SOMETHING TRUE)
      endif()
      
    endforeach(UMFPACK_DEPENDENCY)
    
    if(UMFPACK_OUTPUT MATCHES "${WORD_REGEX}SuiteSparse_")
      find_library(UMFPACK_suitesparseconfig_LIB NAMES suitesparseconfig)
      
      if(NOT UMFPACK_suitesparseconfig_LIB)
        message(FATAL_ERROR "Can not find suitesparseconfig library.")
      endif()
      
      message("!!! Appending ${UMFPACK_suitesparseconfig_LIB} for suitesparseconfig")
      set(UMFPACKLIBS ${UMFPACKLIBS} ${UMFPACK_suitesparseconfig_LIB})
      set(APPENDED_SOMETHING TRUE)
    endif()
    
    message("+++ APPENDED_SOMETHING : ${APPENDED_SOMETHING}")
    
    if(NOT APPENDED_SOMETHING)
      break()
    endif()
    
  endforeach(UMFPACK_TESTS)
  
  if(NOT UMFPACK_COMPILES)
    message(FATAL_ERROR "UMFPACK error")
  else()
    message(STATUS "*** UMFPACK library found")
  endif()
  
endif()

################################################################################
# PETSC

if(PETSC)
  set(PETSC_FOUND NO)
  
  find_package(PkgConfig)
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(PETSC PETSc)
    if(PETSC_FOUND)
      message("+++ PETSC_LDFLAGS : ${PETSC_LDFLAGS}")
      add_definitions(-DPETSC)
    else()
      message(WARNING "*** pkg_check_modules(... PETSc) error : ${PETSC_PKG_ERRORS}")
      endif()
    endif()
  
  if(NOT PETSC_FOUND)
    find_library(PETSC_LDFLAGS NAMES petsc     HINTS ${LOCAL_INSTALL}/lib HINTS "/home/softs/petsc/")
    
    if(PETSC_LDFLAGS)
      add_definitions(-DPETSC)
      message(STATUS "*** PETSC library found")
    else()
      message(WARNING "*** PETSC library not found")
      set(PETSC_LDFLAGS "")
    endif()
    
  endif()
  set(PETSC_LDFLAGS ${PETSC_LDFLAGS} -ldl -lX11)
  message("+++ PETSC_LDFLAGS : ${PETSC_LDFLAGS}")
endif()

################################################################################
# MPI

# cmake may be schizophrenic on MPI variable settings
if(MPI)
  
  # hwloc is necessary for TUGOm
  find_package(PkgConfig)
  if(PKG_CONFIG_FOUND)
    
    pkg_check_modules(HWLOC hwloc)
    if(HWLOC_FOUND)
      message("+++ HWLOC_LDFLAGS : ${HWLOC_LDFLAGS}")
      endif()
    
    endif()
  
  find_package(MPI)
  
  if(MPI_C_COMPILER)
#     set(CMAKE_C_COMPILER ${MPI_CXX_COMPILER})
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  #elseif(MPI_COMPILER) #MPI_COMPILER is mpic++, that fails on SOME versions!!!
  #  set(CMAKE_C_COMPILER ${MPI_COMPILER})
  else()#fail safe
    set(CMAKE_CXX_COMPILER mpicxx)
  endif()
  
  message("+++ CMAKE_CXX_COMPILER : ${CMAKE_CXX_COMPILER}")
  
  if(MPI_C_LINK_FLAGS AND MPI_C_LIBRARIES)
    set(MPI_FLAGS ${MPI_CXX_LINK_FLAGS} ${MPI_CXX_LIBRARIES})
  elseif(MPI_LINK_FLAGS AND MPI_EXTRA_LIBRARY)
    set(MPI_FLAGS ${MPI_LINK_FLAGS} ${MPI_EXTRA_LIBRARY})
  else()#fail safe
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -showme:link OUTPUT_VARIABLE MPI_FLAGS)
  endif()
  
  string(STRIP "${MPI_FLAGS}" MPI_FLAGS)
  
  message("+++ MPI_FLAGS : ${MPI_FLAGS}")
  add_definitions(-DHAVE_MPI)
endif()

################################################################################
# conclusion

#SAVE FLAGS to poc-solvers.LDFLAGS
set(LIBS ${CSPARSELIB} ${HIPSLIBS} ${PASTIXLIBS} ${ARPACKLIB} ${PARPACKLIB} ${MUMPSLIBS} ${SCALAPACKLIB} ${UMFPACKLIBS} ${LAPACKLIB} ${PETSC_LDFLAGS} ${GFORTRANLIB})
string(STRIP "${LIBS}" LIBS)
set(LDFLAGSNAME poc-solvers.LDFLAGS)
set(LDFLAGSPATH ${CMAKE_BINARY_DIR}/${LDFLAGSNAME})
set(LDFLAGSCONTENT ${LIBS} ${MPI_FLAGS} ${HWLOC_LDFLAGS})
execute_process(COMMAND echo ${LDFLAGSCONTENT} OUTPUT_FILE ${LDFLAGSPATH})
message("+++ ${LDFLAGSNAME} : ${LDFLAGSCONTENT}")
#using install to ensure lib is actually created
install(FILES ${LDFLAGSPATH} DESTINATION lib)

#MERCURIAL VERSION
set(REVISION_FILE revision.h)
add_custom_target(${REVISION_FILE} ALL
# creating a temporary file
#  so that the file is only updated if and only if necessary
#   so that targets are remade only when is necessary
  which hg
    && hg -R ${poc-solvers_SOURCE_DIR}/.. parent --template "#define REVISION \"Mercurial revision {rev}:{node|short} of {date|isodate}\"\\n" > ${REVISION_FILE}~
    || touch ${REVISION_FILE}~
  COMMAND
# do not overwrite a non-empty target with an empty temporary file
# and do not touch the target if there is no difference with the temporary file
  test -s ${REVISION_FILE} -a "!" -s ${REVISION_FILE}~
    || diff -q ${REVISION_FILE}~ ${REVISION_FILE}
    || mv -v ${REVISION_FILE}~ ${REVISION_FILE}
    && cp ${REVISION_FILE} ${poc-solvers_SOURCE_DIR}
    || true
  VERBATIM)
include_directories(${CMAKE_BINARY_DIR})

include_directories(${poc-solvers_SOURCE_DIR})

# Only static libraries can be found by the kernel if out of a standard path or LD_LIBRARY_PATH
add_library(poc-solver STATIC initialise.cpp d-utilities.cpp z-utilities.cpp d-factorize.cpp z-factorize.cpp d-solve.cpp z-solve.cpp convert.cpp convertz.cpp read_mat.cpp 
            poc-solvers.cpp mpi-reducers.cpp dynamic.cpp functions.cpp solvers-statistics.cpp matrix.cpp matrix-format.cpp mpi-communications.cpp mpi-utilities.cpp linear-gmres.cpp 
            pastix-interface.cpp hips-interface.cpp mumps-interface.cpp umfpack-interface.cpp cxsparse-interface.cpp lapack-interface.cpp diagonal-interface.cpp gmres-interface.cpp petsc-interface.cpp)

cmake_minimum_required(VERSION 2.6) #because of -lcxsparse
add_executable(testlib testlib.cpp)
target_link_libraries(testlib poc-solver ${TESTLIBS} ${LIBS}) 
add_executable(testlibz testlibz.cpp)
target_link_libraries(testlibz poc-solver ${TESTLIBS} ${LIBS})

if(OPENMP_FOUND)
  set_property(TARGET poc-solver testlib testlibz PROPERTY LINK_FLAGS ${OpenMP_CXX_FLAGS})# -fopenmp
endif()

#INSTALL STUFF
set(CMAKE_INSTALL_PREFIX ".")
install(TARGETS poc-solver testlib testlibz
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)
install(FILES cs.h matrix_io.h poc-solvers.h mpi-reducers.h solvers-functions.h solvers-statistics.h solverlib.h UFconfig.h DESTINATION include)

#MAKE BEHAVIOUR
set(CMAKE_VERBOSE_MAKEFILE TRUE CACHE BOOL "echo all commands to the console during the make" FORCE)

add_custom_target(dist
  ln -s ${poc-solvers_SOURCE_DIR}/.. poc-solvers
  COMMAND tar chf poc-solvers-src.tar.gz poc-solvers/{src/{*.{h,c},CMakeLists.txt},README.{sh,html}}
  COMMAND rm poc-solvers)
