FC = gfortran
LD = gfortran
FFLAGS =-g -O3 -finit-real=zero -mcmodel=medium -fdefault-real-8 -ffpe-trap=invalid,zero,overflow -fopenmp -fbounds-check -fno-align-commons 
#-pedantic
LDFLAGS = `nc-config --flibs` `nc-config --fflags` -fopenmp
OTHERFLAGS = $(LDFLAGS)
#LDFLAGS =-Wl,--no-relax `nc-config --flibs` `nc-config --fflags` -fopenmp
#OTHERFLAGS = -L/opt/local/lib -lnetcdf -lnetcdff -I/opt/local/include


