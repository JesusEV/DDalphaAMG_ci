#
#  This file is part of MUMPS 5.4.0, released
#  on Tue Apr 13 15:26:30 UTC 2021
#
topdir = /home/leemhuis/installs/MUMPS_5.4.0
libdir = $(topdir)/lib

default: c_example


.PHONY: default clean
.SECONDEXPANSION:


include $(topdir)/Makefile.inc

LIBMUMPS_COMMON = $(libdir)/libmumps_common$(PLAT)$(LIBEXT)


LIBDMUMPS = $(libdir)/libdmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)

c_example:	$(LIBDMUMPS) $$@.o
	$(FL) -o $@ $(OPTL) $@.o $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)


.SUFFIXES: .c .F .o
.F.o:
	$(FC) $(OPTF) -I. -I$(topdir)/include -I$(topdir)/src $(INCS) -c $*.F $(OUTF)$*.o
.c.o:
	$(CC) $(OPTC) $(CDEFS) -I. -I$(topdir)/include -I$(topdir)/src $(INCS) -c $*.c $(OUTC)$*.o


$(libdir)/libdmumps$(PLAT)$(LIBEXT):
	@echo 'Error: you should build the library' $@ 'first'
	exit 1

$(LIBMUMPS_COMMON):
	@echo 'Error: you should build the library' $@ 'first'
	exit 1

clean:
	$(RM) *.o c_example 
