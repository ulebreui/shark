#Name of the executable
EXEC=shark

#Compiler
COMP="GNU"

#Number of ghost cells
NGHOST=2

#Number of dust species
NDUST=0
#Geometry (0 == Cartesian, 1 == spherical,2 == 2D-disk)
GEOM=0
GRIDSPACE=0 #0 = linear, 1 = log ! Only for 2D disk
TURB=0

MHD = 0
NDUSTPSCAL = 0 # number of dust passive scalars
GRAVITY    = 0
#Solver 
SOLVER=2 # 0 =lff, 1 = hll, 2 = hllc
#Solver dust
SOLVERDUST=0 #0=Huang & Bai, 1=LLF, 2=HLL
#Solver induction equation
SOLVERB=1 #0=LLF, 1=H&B, 2=HLL
#Number of zones
NX=128 # In X direction
NY=1   # In Y direction
NHBINS=128
# Compiles with open MP
OPENMP=0

#Debug option
DEBUG=0
#Setup used. You can do collapse or 1 zone models. You can implement new setups in the setup folder
SETUP = 
DEFINES = -DNDUST=$(NDUST) -DNX=$(NX) -DNY=$(NY) -DGRIDSPACE=$(GRIDSPACE) -DGEOM=$(GEOM) -DMHD=$(MHD) -DOPENMP=$(OPENMP) -DSOLVER=$(SOLVER) -DSOLVERDUST=$(SOLVERDUST) -DSOLVERB=$(SOLVERB) -DNGHOST=$(NGHOST) -DTURB=$(TURB) -DGRAVITY=$(GRAVITY) -DNDUSTPSCAL=$(NDUSTPSCAL)

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
VPATH  = ../src:../src/commons:../src/grid:../src/turb:../src/solver:../src/setups/$(SETUP):../src/gravity:../src/hydro:../src/chemistry:../src/:../src/dust
#############################################################################
# All objects
MODOBJ = precision.o phys_const.o hydro_commons.o setup_commons.o gravity_commons.o ionisation_commons.o random.o 
ifneq ($(NDUST),0)
MODOBJ+=dust_commons.o 
endif
ifeq ($(TURB),1)
MODOBJ+=turb_commons.o
endif
MODOBJ+= commons.o parameters.o units.o boundary_types.o slope_limiters.o solvers.o allocate_init.o ionisation_utils.o ionisation.o grid.o eos.o setup.o courant.o godunov.o boundaries.o write_output.o read_output.o source_terms.o force_kicks.o fargo_scheme.o
ifneq ($(NDUST),0)
MODOBJ+=dust_init.o dust_drag.o dust_growth.o
endif
OBJS= solve.o read_params.o
ifeq ($(TURB),1)
MODOBJ+= turb_utils.o turb.o
endif
MODOBJ+= gravity_utils.o
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
