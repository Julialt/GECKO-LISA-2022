#---------------------------------------------
# declaration of definitions (macros)

# remember to delete previous *.o files when you change the compiler
#FC = pgf90
FC = gfortran

MODE=PROD
#MODE=DEVEL

ifeq    ($(MODE),DEVEL)
FFLAGS = -g -O0 -fbounds-check -fbacktrace -finit-integer=-100000 -ffpe-trap=zero,overflow,invalid -Waliasing -Wampersand -Wconversion -Wsurprising -Wintrinsics-std -Wintrinsic-shadow -Wline-truncation -Wunused-variable -Wuninitialized
else
FFLAGS = 
endif

#INC = -I.  -I/usr/local/lib
#---------------------------------------------
# VPATH specifies directories to search, in order: VPATH = dir1:dir2:dir3
# so that dependency file.o : file.f is interpreted file.o : dir1/file.f
# vpath (lower case) allows selective specification, e.g.
# vpath %.f ../LIB means files *.f in dir ../LIB

VPATH = ../LIB

INCLUDES = 

EXC = cm

MODULE = keyparameter.mod keyflag.mod minidict.mod \
         tempoci.mod \
         references.mod tempflag.mod database.mod \
         dictstackdb.mod toolbox.mod searching.mod mapping.mod \
         primetool.mod rjtool.mod stdgrbond.mod atomtool.mod brchtool.mod stdtool.mod \
         stdratings.mod ringtool.mod normchem.mod \
         sortstring.mod reactool.mod bensontool.mod tweettool.mod rdkratetool.mod \
         cdtool.mod fragmenttool.mod\
         simpoltool.mod nannoolaltool.mod myrdaltool.mod khydtool.mod gromhetool.mod \
         rdchemprop.mod namingtool.mod switchisomer.mod rxwrttool.mod dictstacktool.mod  \
         spsptool.mod rdmecatool.mod loaddbtool.mod loadc1tool.mod \
         loadchemin.mod masstranstool.mod outtool.mod radchktool.mod dectool.mod \
         roisotool.mod rorxtool.mod rodectool.mod rochem.mod \
         ho_aromtool.mod hotool.mod hoaddtool.mod hochem.mod \
         no3addtool.mod no3abstool.mod no3chem.mod \
         o3addtool.mod poztool.mod criegeetool.mod hotcriegeechem.mod stbcriegeechem.mod panchem.mod o3chem.mod \
         tempo.mod rooohchem.mod ro2tool.mod ro2chem.mod rco3tool.mod rco3chem.mod stuff4hvdiss2.mod \
         dhftool.mod logtool.mod maintool.mod 
         
  
OBJS = rdfixnam.o \
     tempoci.o \
     hvdiss2.o\
     wrt_kc1.o \
     keyparameter.o keyflag.o minidict.o references.o database.o toolbox.o \
     dictstackdb.o dictstacktool.o mapping.o searching.o reactool.o \
     primetool.o rjtool.o stdgrbond.o atomtool.o brchtool.o stdtool.o  stdratings.o \
     ringtool.o normchem.o \
     sortstring.o bensontool.o  tweettool.o rdkratetool.o \
     simpoltool.o nannoolaltool.o myrdaltool.o gromhetool.o khydtool.o rdchemprop.o \
     spsptool.o rdmecatool.o loaddbtool.o loadc1tool.o loadchemin.o \
     rxwrttool.o cdtool.o fragmenttool.o tempflag.o masstranstool.o outtool.o radchktool.o dectool.o\
     roisotool.o rorxtool.o rodectool.o  rochem.o  \
     hoaddtool.o ho_aromtool.o hotool.o hochem.o \
     no3addtool.o no3abstool.o no3chem.o \
     o3addtool.o poztool.o criegeetool.o hotcriegeechem.o stbcriegeechem.o panchem.o o3chem.o \
     ro2tool.o tempo.o rooohchem.o ro2chem.o rco3tool.o rco3chem.o stuff4hvdiss2.o switchisomer.o \
     namingtool.o dhftool.o logtool.o maintool.o main.o

#---------------------------------------------
# relationship (dependency) between executable (target) & object files,
#       command to compile executable
#       $@ means 'current target'
$(EXC) : $(MODULE) $(OBJS) 
	$(FC) -fdefault-real-8 -mcmodel=large $(FFLAGS) $(OBJS) -o $@ 
	echo $(LIB)
	mv $@ ../RUN

#---------------------------------------------
# relationship (dependency) between object files (target) and dependent files
$(OBJS) : $(INCLUDES)

# make the modules 
%.mod : %.f90
	$(FC) -fdefault-real-8 -mcmodel=large -c $(FFLAGS) -s $<
	touch $*.mod $*.o

#---------------------------------------------
# test f90
%.o : %.f90 
	$(FC) $(FFLAGS) -fdefault-real-8 -mcmodel=large -c $(INC) $<
#	gfortran -mcmodel=medium -c -s $<


#---------------------------------------------
# suffixes rule: instructs compiler to find files with any undernoted suffix
.SUFFIXES : .f .h .o

#---------------------------------------------
# inference rule: compile all .f files into .o files
#       $* means 'current target without its extension'
#       $< means 'first dependent file, including directory names'
#       $^ means 'all dependents, including directory names' 
.f.o :  
	$(FC) $(FFLAGS) -fdefault-real-8 -mcmodel=large -c $(INC) $<

#---------------------------------------------
# relationship (dependency) between object (target) and fortran files
# are not essential because of the .f.o inference rule!

# main.o   :   main.f  
#hvdiss2.o :   hvdiss2.f    

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod
