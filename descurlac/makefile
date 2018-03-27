CC   = gcc
PREFIX = $(shell pwd)
SRCDIR = $(PREFIX)/src
OBJDIR = $(PREFIX)/obj
MODDIR = $(PREFIX)/mod
LIBDIR = $(PREFIX)/lib
SRC  := $(wildcard $(SRCDIR)/*.f90)
OBJ2  := $(patsubst $(SRCDIR)%.f90, $(OBJDIR)%.o, $(SRC))
OBJ = $(OBJDIR)/mod_konst.o $(OBJDIR)/mod_curve3d.o $(OBJ2)
EXE = $(EXEDIR)/DESCUR_SUN
SOBJ = $(LIBDIR)/descur_idl.so
MODEXT = mod
DSCOBJ = $(OBJDIR)/descur_idl.o
ARC = $(LIBDIR)/libdescur.a
MFLAGS = -M
MODFLAGS = $(MFLAGS)$(MODDIR)

system = $(shell uname -s)
hardware = $(shell uname -m)
ifeq ($(system),Linux)
	AR   = /usr/bin/ar r
	FC90 = gfortran
#	FFLAGS =  -fpic -c
#	CFLAGS =  -fpic -c
#	LDFLAGS = -fpic 
	FFLAGS =  -fpic -c
	CFLAGS =  -fpic -c
	LDFLAGS = -fpic 
endif
# ifeq ($(system),Linux)
# 	AR   = /usr/bin/ar r
# 	FC90 = ifort
# 	FFLAGS =  -fpic -c
# 	CFLAGS =  -fpic -c
# 	LDFLAGS = -fpic -fpp
# endif


# Sequence is important for .so library
OBJECTS = mod_konst mod_curve3d amplitud evolve fixarray fftrans funct  \
          getangle inguess restart scrunch curve3d readin order
OBJ_SET := $(patsubst %, $(OBJDIR)/%.o, $(OBJECTS) )

# $(SOBJ): $(DSCOBJ)
# 	$(FC90) $(LDFLAGS) -G -z defs -check all -o $(SOBJ) $(DSCOBJ) $(OBJ_SET) $(MATHLIBS)
# 	$(AR) $(ARC) $(OBJDIR)/*.o
# 	@mv *.$(MODEXT) $(MODDIR)/
$(SOBJ): $(DSCOBJ)
	$(FC90) $(LDFLAGS) -g -z defs -shared -o $(SOBJ) $(DSCOBJ) $(OBJ_SET) $(MATHLIBS)
	$(AR) $(ARC) $(OBJDIR)/*.o
	@mv *.$(MODEXT) $(MODDIR)/

# ===============
# Implicit rules:
# ===============

$(DSCOBJ):  $(OBJ) $(SRCDIR)/descur_idl.c
	$(CC) $(CFLAGS) $(SRCDIR)/descur_idl.c -o $(DSCOBJ)

$(OBJDIR)/mod_konst.o: $(SRCDIR)/mod_konst.F90
	$(FC90) $(FFLAGS) $(SRCDIR)/mod_konst.F90 -o $(OBJDIR)/mod_konst.o

$(OBJDIR)/mod_curve3d.o: $(SRCDIR)/mod_curve3d.F90
	$(FC90) $(FFLAGS) $(SRCDIR)/mod_curve3d.F90 -o $(OBJDIR)/mod_curve3d.o

$(OBJDIR)%.o: $(SRCDIR)%.f90
	$(FC90) $(FFLAGS) $< -o $@
