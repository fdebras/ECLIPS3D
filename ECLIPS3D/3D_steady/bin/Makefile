#=====================================================
# Run directory file
#---------------------------------------------------
RUN = ../run/
#
#---------------------------------------------------
#=====================================================


#############################################################################
#---------------------------------------------------
F90 = mpif90
FLAGS = 
CPPFLAGS = 
#FLAGS    =  -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include#-check_mpi
#FLAGS    = -ffast-math -ffree-line-length-0 -fdefault-real-8  \
#    -fbounds-check -ffpe-trap=invalid,zero,overflow,underflow \
#    -g -fbacktrace -Warray-bounds
# Libraries
LIBS     = /usr/lib/libblacs-openmpi.so.1 /usr/lib/libscalapack-openmpi.so.1.8.0 \
    /usr/lib/libblacsF77init-openmpi.so.1
# -lblas -lcurl 
#---------------------------------------------------

#---------------------------------------------------
# directories and vpath for source files
#---------------------------------------------------
BIN = ../bin/
SRC = ../src/
# VPATH  = $(CASDIR):$(SRC_3D):$(SRC_2D_axi):$(SRC_2D_shallow)
#---------------------------------------------------
# Source objects
#---------------------------------------------------
SRCOBJ = mod_data.o mod_init_para.o \
	mod_init_matrix.o mod_3D_fill_matrix.o writing_test.o mod_solver.o

###########################################################################
#############################################################################	

3D_steady :		$(SRCOBJ) standing_solver_3D.o
				$(F90) $(FLAGS) $(SRCOBJ) standing_solver_3D.o -o \
				$(RUN)ECLIPS3D_steady.exe $(LIBS)  

standing_solver_3D.o:		$(SRC)standing_solver_3D.F90 Makefile
		$(F90) $(CPPFLAGS) $(FLAGS) -c $(SRC)standing_solver_3D.F90
				
%.o:			$(SRC)%.F90
		$(F90) $(CPPFLAGS) $(FLAGS) -c $< -o $@


#############################################################################


#############################################################################
clean :
		rm  *.o *.mod

#############################################################################

