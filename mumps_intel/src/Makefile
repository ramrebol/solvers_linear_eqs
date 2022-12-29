# MUMPS was compiled using ifort and MKL
# Created            : 15/07/2017
# Last modificacion  : 28/12/2022

include template

.DEFAULT:
	-touch $@
all:run
decimal.o: ./decimal.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c ./decimal.f90
types.o: ./decimal.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c ./types.f90
mumps_system.o: ./mumps_system.f90 decimal.o types.o
	$(FC) $(OPTL) -L$(LIBMUMPS)/libseq -I. -I$(IMUMPS) -c ./mumps_system.f90
system.o: ./system.f90 decimal.o types.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c ./system.f90
main.o: main.f90 decimal.o types.o mumps_system.o system.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c ./main.f90

SRC = ./main.f90 ./decimal.f90 ./types.f90 ./mumps_system.f90 ./system.f90

OBJ = main.o decimal.o types.o mumps_system.o system.o

clean: neat
	-rm -f $(OBJ) ../bin/run *.mod
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
run: $(OBJ)
	$(LD) -o ../bin/run $(OPTL) $(OBJ) $(LIBSMKL) $(LIBSMUMPS)
