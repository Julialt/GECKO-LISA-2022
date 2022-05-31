FC = gfortran
LD = gfortran
FFLAGS = -O3 -ffast-math -march=native -fdefault-real-8 -finit-real=nan -Wall -std=f2008ts -fall-intrinsics `nc-config --fflags`
LDFLAGS = ./geckolib.a `nc-config --flibs`

