#=====================================================
# Run directory file
#---------------------------------------------------
RUN = ../run/
#
#---------------------------------------------------
#=====================================================


#############################################################################
#---------------------------------------------------
# Compiler and flags
F90 = mpif90
FLAGS =
CPPFLAGS =
#FLAGS = -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include#-check_mpi
#FLAGS = -ffast-math -ffree-line-length-0 -fdefault-real-8  \
#        -fbounds-check -ffpe-trap=invalid,zero,overflow,underflow \
#        -g -fbacktrace -Warray-bounds
#---------------------------------------------------

#---------------------------------------------------
# Libraries
LIBS = /usr/lib/libblacs-openmpi.so.1 /usr/lib/libscalapack-openmpi.so.1.8.0 \
       /usr/lib/libblacsF77init-openmpi.so.1
#---------------------------------------------------

#---------------------------------------------------
# Directories and vpath for source files
BIN = ../bin/
SRC = ../src/
VPATH = $(SRC)
#---------------------------------------------------

#---------------------------------------------------
# Source files
#---------------------------------------------------

SRCFILES = mod_3D_fill_matrix.F90 mod_data.F90 mod_eigenvalues.F90 \
           mod_init_matrix.F90 mod_init_para.F90 mod_writing.F90 \
           normal_modes_3D.F90

OBJFILES := $(SRCFILES:.F90=.o)

TARGET = $(RUN)ECLIPS3D.exe
###########################################################################
all: $(TARGET)

%.o:%.F90
	$(F90) $(CPPFLAS) $(FLAGS) -c $< -o $@

clean:
	@rm *.o *.mod

$(TARGET) : $(OBJFILES)
	$(F90) $(FLAGS) $(OBJFILES) -o $(TARGET) $(LIBS)

mod_3D_fill_matrix.o: mod_3D_fill_matrix.F90 mod_data.o mod_init_matrix.o \
  mod_init_para.o
mod_data.o: mod_data.F90
mod_eigenvalues.o: mod_eigenvalues.F90 mod_data.o mod_3D_fill_matrix.o \
  mod_writing.o
mod_init_matrix.o: mod_init_matrix.F90 mod_data.o mod_init_para.o
mod_init_para.o: mod_init_para.F90
mod_writing.o: mod_writing.F90
normal_modes_3D.o: normal_modes_3D.F90 mod_3D_fill_matrix.o mod_data.o \
  mod_eigenvalues.o mod_init_matrix.o mod_init_para.o
#############################################################################
