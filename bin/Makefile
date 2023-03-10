#Name of the executable
EXEC=shark

#Compiler
COMP="GNU"

#Number of ghost cells
NGHOST=2

#Number of dust species
NDUST=0
#Use of energy equation (not coded yet)
ENERGY=0
#Use of full MHD (not coded yet) / Otherwise approx of Lebreuilly et al., 2023
MHD=0
#Gravity (assumes spherical geometry)
GRAVITY=0
#Spherical geometry (0 == Cartesian)
SPHERE=0
#Number of zones
NHBINS=128
# Compiles with open MP
OPENMP=0

#Debug option
DEBUG=0
#Setup used. You can do collapse or 1 zone models. You can implement new setups in the setup folder
SETUP = collapse_spher
DEFINES = -DNDUST=$(NDUST) -DNHBINS=$(NHBINS) -DENERGY=$(ENERGY) -DMHD=$(MHD) -DGRAVITY=$(GRAVITY) -DSPHERE=$(SPHERE) -DOPENMP=$(OPENMP) -DNGHOST=$(NGHOST)

#############################################################################
# Fortran compiler options and directives
#############################################################################
# GFORTRAN
ifeq ($(COMP),"GNU")
F90 = gfortran 
ifeq ($(OPENMP),1)
F90 = gfortran -fopenmp 
endif
   FFLAGS +=  -x f95-cpp-input -O3 -frecord-marker=4 -ffree-line-length-none $(DEFINES)
ifeq ($(DEBUG),1)
   FFLAGS += -fbounds-check -fbacktrace -Wuninitialized -Wall
   FFLAGS += -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan
   FFLAGS += -finit-real=nan
endif
endif
#INTEL
ifeq ($(COMP),"INTEL")
F90 = ifort
ifeq ($(OPENMP),1)
F90 = ifort -qopenmp
endif
FFLAGS = -cpp $(DEFINES) -DWITHOUTMPI
FFLAGS += -O3 -fp-model source 
endif

#############################################################################m
# Sources directories are searched in this exact order
BINDIR = .
VPATH =../src:../src/commons:../src/solver:../src/setups/$(SETUP):../src/hydro:../src/gravity:../src/dust:../src/mhd
#############################################################################
# All objects
MODOBJ = precision.o phys_const.o hydro_commons.o mhd_commons.o setup_commons.o
ifneq ($(NDUST),0)
MODOBJ+=dust_commons.o 
endif
ifneq ($(GRAVITY),0)
MODOBJ+=gravity_commons.o
endif
MODOBJ+= commons.o parameters.o units.o allocate_init.o eos.o B_field.o mhd_utils.o setup.o courant.o godunov.o boundaries.o write_output.o read_output.o source_terms.o
ifneq ($(GRAVITY),0)
MODOBJ+=gravity_utils.o
endif
ifneq ($(NDUST),0)
MODOBJ+=dust_init.o dust.o dust_evol.o charge.o
endif
OBJS= solve.o read_params.o

#############################################################################
all:
	for file in $(MODOBJ); do unset MAKELEVEL ; $(MAKE) $$file; done
	$(MAKE) shark
shark: $(MODOBJ) $(OBJS) shark.o
	   $(F90) $(MODOBJ) $(OBJS) shark.o -o ./$(EXEC)

%.o:%.F
	$(F90) $(FFLAGS) -c $^ -o $@ $(LIBS_OBJ)
%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@ $(LIBS_OBJ)
%.o:%.f95
	$(F90) $(FFLAGS) -c $^ -o $@ $(LIBS_OBJ)
FORCE:
#############################################################################
clean :
	rm -f *.o *.$(MOD) *__genmod.f90 *mod shark
#############################################################################
