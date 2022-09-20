
VPATH=../../src


STANDART_SOFTS=/sirocco/ocean/softs
LOCAL_SOFTS=/home2/ocean/softs/gnu-20130108

#STANDART_SOFTS=/home/softs
#LOCAL_SOFTS=/home/softs

#==========================================================
# DEFINE OPTIONAL LIBRARIES MUMPS UMFPACK LAPACK
#  -------- MUMPS --------------------
MUMPS=YES
MUMPSDIR   = /home/ocean/softs/gnu/Solver-dev/MUMPS_5.0.1/
#$(STANDART_SOFTS)/MUMPS_4.9.2/parallel
# --------- UMPFPACK ------------------
UMFPACK=YES
UMFPACKDIR   = $(STANDART_SOFTS)/SuiteSparse/UMFPACK
AMDDIR       = $(STANDART_SOFTS)/SuiteSparse/AMD
UMSCONFIGDIR = $(STANDART_SOFTS)/SuiteSparse/UFconfig/
#  -------- HYPRE --------------------
HYPRE=NO
HYPREDIR   = $(STANDART_SOFTS)/hypre-2.7.0b/
# ---------- HIPS --------------------
HIPS=YES
HIPSDIR    = /home/ocean/softs/gnu/Solver-dev/zdhips
#$(STANDART_SOFTS)/hips-complex/LIB/
SCOTCHDIR  = /home/ocean/softs/gnu/Solver-dev/scotch_6.0.4/lib
#$(STANDART_SOFTS)/scotch_6.0.0/lib/
# --------- LAPACK --------------------
LAPACK=NO
LAPACKLIB=
# LAPACK is include in blas librarie with acml
#  END OPTIONAL SECTION 
# PROF
#PROF = -g -pg
# ---------- PASTIX --------------------
PASTIX=YES
#  HIPSDIR    = 
#  SCOTCHDIR  = 
## --------- LAPACK --------------------
# --------- SCALAPACK --------------------
SCALAPACKLIB= $(STANDART_SOFTS)/scalapack/lib/
#==========================================================

# BLAS (fortran one)
LIBBLAS = -lblas ${LOCAL_SOFTS}/scalapack/lib/blacs.a  ${LOCAL_SOFTS}/scalapack/lib/blacsF77.a ${LOCAL_SOFTS}/scalapack/lib/blacs.a 

LIBOTHERS = -lgfortran -lpthread -lm -lc

ifeq ($(MUMPS),YES)
	INC1M = -I$(MUMPSDIR)/include 
	LOCAL_CFLAGSM = -DMUMPS -DHAVE_MPI
	OBJS1M = $(MUMPSDIR)/lib/libdmumps.b \
	         $(MUMPSDIR)/lib/libpord.b
        MUMPSMLIB=$(MUMPSDIR)/lib/libdmumps.a
        MUMPSCLIB=$(MUMPSDIR)/lib/libmumps_common.a
        MUMPSPLIB=$(MUMPSDIR)/PORD/lib/libpord.a
#        MUMPSSLIB=$(MUMPSDIR)/libseq/libmpiseq.a 
        MUMPSLIBT = $(MUMPSDIR)/lib/libdmumps.a $(MUMPSDIR)/lib/libmumps_common.a -L/home/ocean/softs/gnu/Solver-dev/metis/lib -lmetis -L$(MUMPSDIR)/PORD/lib/ -lpord -L/home/ocean/softs/gnu/Solver-dev/scotch_6.0.4/lib -lptesmumps -lptscotch -lptscotcherr -lscotch -L/home/ocean/softs/gnu/Solver-dev/scalapack-2.0.2/ -lscalapack   -lblas -llapack -lpthread

endif

ifeq ($(UMFPACK),YES)
	INC1U = -I$(UMFPACKDIR)/Include -I$(AMDDIR)/Include	-I$(UMSCONFIGDIR)	
	LOCAL_CFLAGSU = -DUMFPACK
	OBJS1U = $(UMFPACKDIR)/Lib/libumfpack.b $(AMDDIR)/Lib/libamd.b
	UMFLIB=$(UMFPACKDIR)/Lib/libumfpack.a
	AMDLIB=$(AMDDIR)/Lib/libamd.a
endif

ifeq ($(HYPRE),YES)
        INC1H = -I $(HYPREDIR)/include
        LOCAL_CFLAGSH = -DHYPRE -DHAVE_MPI
        HYPRELIB = $(HYPREDIR)/lib/libHYPRE.a
endif

ifeq ($(HIPS),YES)
        INC1H = -I$(HIPSDIR)/DLIB -I$(HIPSDIR)/ZLIB
        LOCAL_CFLAGSH = -DHIPS -DHAVE_MPI 
        HIPSLIB1 = -L$(HIPSDIR)/DLIB -ldhips  -ldio -lm -ldspkit  -L$(SCOTCHDIR) -lscotch -lscotcherr
        HIPSLIB = $(HIPSLIB1) -L$(HIPSDIR)/ZLIB -lzhips  -lzio -lm -lzspkit  -L$(SCOTCHDIR) -lscotch -lscotcherr
# 	HIPSLIB1=$(HIPSDIR)/libhips.a
# 	HIPSLIB2=$(HIPSDIR)/libio.a
# 	HIPSLIB3=$(HIPSDIR)/libspkit.a
# 	HIPSLIB0=$(SCOTCHDIR)/libscotchmetis.a
# 	HIPSLIB4=$(SCOTCHDIR)/libscotch.a
# 	HIPSLIB5=$(SCOTCHDIR)/libscotcherr.a
endif

ifeq ($(PASTIX),YES)
        INC1P = -I/home/ocean/softs/gnu/Solver-dev/pastix/include
        LOCAL_CFLAGSP = -DDISTRIBUTED -DWITH_SCOTCH -DPASTIX -DHAVE_MPI 
#-DINTSIZE32  -DPREC_DOUBLE -DTYPE_COMPLEX -DFORCE_COMPLEX -DTYPE_COMPLEX -DPREC_DOUBLE  -DCUDA_SM_VERSION=20 -DMEMORY_USAGE -DSTATS_SOPALIN -DPASTIX_DYNSCHED -I/sirocco/ocean/softs/pastix_release_4030_parallel/include
#        PASTIXLIB1 = -L/home2/ocean/softs/Gnu/pastix_release_4030/src/../install -lpastix -lmatrix_driver -lpastix -lgfortran -lm -lrt -L/sirocco/ocean/softs/metis-4.0/ -lmetis -L/home2/ocean/softs/Gnu/scotch_6.0.0//lib -lptscotch -lscotch -lptscotcherrexit -lpthread
	 PASTIXLIB1 = -L/home/ocean/softs/gnu/Solver-dev/pastix/lib -lpastix_murge -lpastix -lmatrix_driver -lgfortran -lm -lrt -L/home/ocean/softs/gnu/Solver-dev/scotch_6.0.4/lib -lptscotch -lscotch -lptscotcherrexit -L/home/ocean/softs/gnu/lib64 -lhwloc -lpthread
endif
ifeq ($(LAPACK),YES)
	LOCAL_CFLAGSL = -DLAPACK
endif

CC=mpicc

INC1 = -I../include/ $(INC1U) $(INC1H) $(INC1P) $(INC1M) $(INC1H)

LOCAL_CFLAGS = $(LOCAL_CFLAGSM) $(LOCAL_CFLAGSU)  $(LOCAL_CFLAGSH) $(LOCAL_CFLAGSL) $(LOCAL_CFLAGSP) $(PROF) -DDEBUG -DVERBOSE -g -fopenmp

OBJS1 = $(OBJS1M) $(OBJS1U) $(OBJS1L)

OBJS =  solver.o convert.o convertz.o read_mat.o

FL = gfortran

default: lib testlib testlibz

all: lib testlib testlibz


.SUFFIXES: .F .o
.F.o:
	$(FC) $(OPTF) $(INC) -I. -I../include -g -c -O3 $*.F

.SUFFIXES: .c .o 
.c.o:
#	$(CC) $(LOCAL_CFLAGS)  $(OPTC) $(INC) $(INC1) -I.  -c -g $<
	$(CC) $(LOCAL_CFLAGS)  $(OPTC) $(INC) $(INC1) -I. -g -fbounds-check -c -O3 $<

.SUFFIXES: .a .b
.a.b:
	ar x $<

clean:
	$(MAKE) -C ../../CXSparse clean
	$(RM) *.o *.a  
	$(RM) ../lib/*.a  
	$(RM) ../include/*  

purge : 
	$(MAKE) -C ../../CXSparse purge
	$(RM) *.o libsolver.a


testlib:  testlib.o ../lib/libpoc-solver.a
	$(CC) -g -O3 -o testlib $(LOCAL_CFLAGS)  $(OPTL) testlib.o -lmpi_f77 -L$(SCALAPACKLIB) ../lib/libpoc-solver.a ../lib/libpoc-solver-domestic.a -lscalapack \
	$(LIBBLAS) $(LIBOTHERS)  $(UMFLIB) $(AMDLIB) $(MUMPSMLIB) $(HIPSLIB) $(PASTIXLIB1) $(MUMPSLIBT)
#	$(SCALAPACKLIB)/blacsC.a $(SCALAPACKLIB)/blacsF77.a $(SCALAPACKLIB)/blacs.a $(SCALAPACKLIB)/blacsF77.a $(SCALAPACKLIB)/blacs.a $(LIBBLAS) $(HYPRELIB) $(HIPSLIB) $(LIBOTHERS) 

testlibz:  testlibz.o ../lib/libpoc-solver.a
	$(CC) -g -O3 -o testlibz $(LOCAL_CFLAGS)  $(OPTL) testlibz.o -L$(SCALAPACKLIB) ../lib/libpoc-solver.a ../lib/libpoc-solver-domestic.a -lmpi_f77 -lmpi -lscalapack \
	$(LIBBLAS)  $(LIBOTHERS) $(UMFLIB) $(AMDLIB) $(MUMPSMLIB) $(HIPSLIB) $(PASTIXLIB1) $(MUMPSLIBT)
#	$(SCALAPACKLIB)/blacsC.a $(SCALAPACKLIB)/blacsF77.a $(SCALAPACKLIB)/blacs.a $(SCALAPACKLIB)/blacsF77.a $(SCALAPACKLIB)/blacs.a $(LIBBLAS)  $(HYPRELIB) $(HIPSLIB)  $(LIBOTHERS) 


lib :   init $(OBJS) $(OBJS1)
	$(MAKE) -C ../../CXSparse library
	ar x ../../CXSparse/Lib/libcxsparse.a
	ar cru libpoc-solver-domestic.a *.o
	mkdir -p ../lib/
	cp libpoc-solver-domestic.a ../lib/
	cp $(VPATH)/*.h ../include/
	echo $(UMFLIB) $(AMDLIB) $(MUMPSMLIB) $(HIPSLIB) $(PASTIXLIB1) $(MUMPSLIBT) > ../lib/poc-solvers.LDFLAGS

ifeq ($(MUMPS),YES)
# 	ar x $(MUMPSMLIB)
# 	ar x $(MUMPSCLIB)
# 	ar x $(MUMPSPLIB)
# #	ar x $(MUMPSSLIB)
endif
ifeq ($(UMFPACK),YES)
# 	ar x $(UMFLIB)
# 	ar x $(AMDLIB)
# 	ar x /usr/lib64/libblas.a
endif
ifeq ($(HIPS),YES)
# 	ar x $(HIPSLIB1)
# 	ar x $(HIPSLIB2)
# 	#ar x $(HIPSLIB3)
# 	ar x $(HIPSLIB0)
# 	ar x $(HIPSLIB4)
# 	ar x $(HIPSLIB5)
endif
	ar cru libpoc-solver.a *.o
	cp libpoc-solver.a ../lib/

init:
	mkdir -p ../include/
ifeq ($(MUMPS),YES)
# 	echo $(MUMPSDIR)
# 	cp $(MUMPSDIR)/include/*.h ../include/
# #	cp $(MUMPSDIR)/libseq/*.h ../include/
# 	cp $(MUMPSDIR)/PORD/include/*.h ../include/
endif
ifeq ($(UMFPACK),YES)
# 	cp $(UMFPACKDIR)/../UFconfig/UFconfig.h ../include/
# 	cp $(UMFPACKDIR)/Include/*.h ../include/
# 	cp $(AMDDIR)/Include/*.h ../include
endif
ifeq ($(HIPS),YES)
# 	cp $(HIPSDIR)/*.h ../include/
endif


