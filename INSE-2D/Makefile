# Directories
HOMEDIR = $(shell pwd | sed -e 's/\/src.*//')
LIBDIR  = $(HOMEDIR)
MODDIR  = $(HOMEDIR)
OBJDIR  = $(HOMEDIR)
BINDIR  = $(HOMEDIR)
VPATH   = $(LIBDIR) $(BINDIR) $(OBJDIR)
OPTDIR  = /WORK2/pku_yyg/mengzy/opt

# Compiler and archiver
CC  = /usr/local/mpi3-dynamic/bin/mpicc
CXX = /usr/local/mpi3-dynamic/bin/mpicxx
F90 = /usr/local/mpi3-dynamic/bin/mpif90
F77 = /usr/local/mpi3-dynamic/bin/mpif90
LD  = /usr/local/mpi3-dynamic/bin/mpif90
AR  = ar rcv
RL  = ranlib

# Compiler flags
CFLAGS   =
F90FLAGS = 
LDFLAGS  = 
INCFLAGS = 
MODFLAGS = 
DBGFLAGS = -g -CA -CB -CS -CV -traceback -debug all -ftrapuv -check all -WB -warn all
OPTFLAGS = -O3
#OPTFLAGS = -O3 -xCORE-AVX-I -ip

FFTW_DIR = $(OPTDIR)/fftw3
FFTW_INC = -I$(FFTW_DIR)/include
FFTW_LIB = -L$(FFTW_DIR)/lib -lfftw3

# Installation script
INSTDIR = $(HOME)/bin
INSTSCPT = cp $(BINDIR)/* $(INSTDIR)/.

F90FILES = DNS.f90 subroutine.f90
# BINFILE  = DNS
OFILES   = $(F90FILES:.f90=.o)
# MODFILES = $(F90FILES:.f90=.mod)

FLGAS = -O3 -CB
#FLAGS = -O3 -fopenmp -CB #-fopenmp -debug all -traceback -ftrapuv -WB -check all#-O3 -openmp -CB #-debug #-g -CA -CB -CS -CV -traceback -debug all -ftrapuv -check all -WB -warn all
#FLAGS = -debug all -traceback -check all 

.SUFFIXES: .f90 .o

comu: $(OFILES)
	cd $(OBJDIR); $(LD) $(FLAGS) -o DNS ${OFILES} ${FFTW_LIB}

.f90.o: 
	$(F90) $(FLAGS) $(INCFLAGS) $(FFTW_INC) -c $*.f90 -o $(OBJDIR)/$*.o $(MODFLAGS)

clean: 
	cd $(OBJDIR); rm -f $(OFILES)
	cd $(MODDIR); rm -f $(MODFILES)
	cd $(BINDIR); rm -f $(BINFILE)

