# Makefile created by mkmf $Id: mkmf,v 18.0 2010/03/02 23:26:08 fms Exp $ 



include template.mk


.DEFAULT:
	-echo $@ does not exist.
all: inca
chkaux.o: ./chkaux.f
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./chkaux.f
chkreac.o: ./chkreac.f
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./chkreac.f
chkspec.o: ./chkspec.f sorting_module.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./chkspec.f
getaux.o: ./getaux.f
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./getaux.f
getreac.o: ./getreac.f akparameter_module.o sorting_module.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./getreac.f
getspec.o: ./getspec.f
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./getspec.f
inca.o: ./inca.f akparameter_module.o sorting_module.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./inca.f
ncutil.o: ./ncutil.f
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./ncutil.f
akparameter_module.o: ./akparameter_module.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./akparameter_module.f90
sorting_module.o: sorting_module.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./sorting_module.f90
keyflag.o: keyflag.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./keyflag.f90
util.o: ./util.f akparameter_module.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./util.f
wrtlinkbin.o: ./wrtlinkbin.f akparameter_module.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./wrtlinkbin.f
wrtlinkncdf.o: ./wrtlinkncdf.f akparameter_module.o sorting_module.o keyflag.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./wrtlinkncdf.f
SRC = ./ncutil.f ./akparameter_module.f90 ./getreac.f ./getspec.f ./wrtlinkbin.f ./inca.f ./getaux.f ./chkspec.f ./sorting_module.f90 ./keyflag.f90 ./chkaux.f ./util.f ./wrtlinkncdf.f ./chkreac.f
OBJ = ncutil.o getreac.o getspec.o wrtlinkbin.o inca.o getaux.o chkspec.o sorting_module.o keyflag.o chkaux.o util.o wrtlinkncdf.o chkreac.o
clean: neat
	-rm -f .inca.cppdefs $(OBJ) inca
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
inca: $(OBJ) 
	$(LD) $(OBJ) -o inca  $(LDFLAGS)
