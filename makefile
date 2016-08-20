SRC_C_ALL:= $(wildcard src_c/*.cpp)
SRC_C_ALL+= $(wildcard src_c/MeshLogic/*.cpp)
SRC_C_ALL+= $(wildcard src_c/Problem/*.cpp)
SRC_C_ALL+= $(wildcard src_c/Solver/*.cpp)

#first we must compile def files and then the source files
SRC_F_ALL:= $(wildcard src_f/def_*.*f90)
SRC_F_ALL+= $(filter-out $(wildcard src_f/def_*.*f90),$(wildcard src_f/*.f90))

OBJECT_ALL:= ./*.o


BINARY:= MeshAndSolver

CC  = g++
FORT = gfortran
LD  = g++

CFLAGS =
FFLAGS = -c -g -ffixed-line-length-none -ffree-line-length-none
CPPFLAGS = -std=c++11 -c -g
LDFLAGS = -lgfortran



all:
	$(FORT) $(FFLAGS) $(SRC_F_ALL)
	$(CC) $(CPPFLAGS) $(SRC_C_ALL) 
	$(LD) -o $(BINARY) $(OBJECT_ALL) $(LDFLAGS)

clean:
	rm ./*.o
	rm ./*.mod
	rm ./$(BINARY)

