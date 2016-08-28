CC=gcc 
#CC=gcc -fopenmp
CCMP=mpicc
#CCMP=/usr/lib64/openmpi/bin/mpicc -fopenmp
#LFLAGS -pg
#LFLAGS= -Wall -pg -fopenmp -lm -L/usr/local/lib
#LFLAGS= -Wall -pg  -lm -L/usr/local/lib
LFLAGS= -lm -Wall -L/opt/openmpi-1.10.1/lib

#LGSLFLAGS= -lgsl -lgslcblas 
CFLAGS= -Wall -I/opt/local/include -Iinclude
OPTFLAG= -O2

SRCDIR=src
OBJDIR=obj
INCL = -I../include
INCMPIFLAG=-D ACTIVE_MPI
INCSMPFLAG=-D STEMPERING

BINDIR=bin


ifeq ($(version),MPI)

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CCMP) $(INCMPIFLAG) $(OPTFLAG) $(CFLAGS) -c $< $(INCL) -o $@

montegrappa_MPI: $(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/mc.o $(OBJDIR)/local_move.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/potential.o $(OBJDIR)/optimizepot.o $(OBJDIR)/montegrappa.o $(OBJDIR)/MPIfunc.o
	$(CCMP) $(INCMPIFLAG) $(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/mc.o $(OBJDIR)/local_move.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/potential.o $(OBJDIR)/optimizepot.o $(OBJDIR)/montegrappa.o $(OBJDIR)/MPIfunc.o -o $(BINDIR)/montegrappa_mpi $(LFLAGS) 

$(OBJDIR)/MPIfunc.o: $(SRCDIR)/MPIfunc.c
	$(CCMP) $(INCMPIFLAG) $(CFLAGS) -c $< $(INCL) -o $@

else ifeq ($(version),STEMPERING) 

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(INCSMPFLAG) $(OPTFLAG) $(CFLAGS) -c $< $(INCL) -o $@


montegrappa_stemp:	$(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/mc.o $(OBJDIR)/local_move.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/potential.o $(OBJDIR)/optimizepot.o $(OBJDIR)/montegrappa.o $(OBJDIR)/stempering.o $(OBJDIR)/memory1.o $(OBJDIR)/adjust_st.o $(OBJDIR)/do_mhistogram.o
	$(CC) $(INCSMPFLAG) $(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/mc.o  $(OBJDIR)/local_move.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/potential.o $(OBJDIR)/optimizepot.o $(OBJDIR)/montegrappa.o $(OBJDIR)/stempering.o $(OBJDIR)/memory1.o $(OBJDIR)/adjust_st.o $(OBJDIR)/do_mhistogram.o  -o $(BINDIR)/montegrappa $(LFLAGS) $(LGSLFLAGS)

else

$(OBJDIR)/%.o: $(SRCDIR)/%.c 
	$(CC) $(OPTFLAG) $(CFLAGS) -c $< $(INCL) -o $@

montegrappa:  $(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/mc.o $(OBJDIR)/local_move.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/potential.o $(OBJDIR)/optimizepot.o $(OBJDIR)/montegrappa.o 
	$(CC) $(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/mc.o  $(OBJDIR)/local_move.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/potential.o $(OBJDIR)/optimizepot.o $(OBJDIR)/montegrappa.o -o $(BINDIR)/montegrappa $(LFLAGS)

endif



grappino:    $(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/potential.o $(OBJDIR)/energy.o $(OBJDIR)/pdb.o $(OBJDIR)/rotamers.o $(OBJDIR)/grappino.o $(OBJDIR)/optimizepot.o
	$(CC) $(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/potential.o $(OBJDIR)/energy.o $(OBJDIR)/pdb.o $(OBJDIR)/rotamers.o $(OBJDIR)/optimizepot.o $(OBJDIR)/grappino.o -o $(BINDIR)/grappino $(LFLAGS)

mhistogram:
	cd src/mhistogram; make

mgp2pdb:     $(OBJDIR)/mgp2pdb.o $(OBJDIR)/io.o $(OBJDIR)/memory.o $(OBJDIR)/geometry.o $(OBJDIR)/misc.o
	$(CC) $(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/mgp2pdb.o -o $(BINDIR)/mgp2pdb $(LFLAGS)

clean:
	rm -f $(OBJDIR)/*.o src/mhistogram/*.o  $(BINDIR)/* 

cleanobj:
	rm -f $(OBJDIR)/*.o 
	rm -f src/mhistogram/*.o



all:
	make cleanobj;
	make version=STEMPERING;
	make cleanobj;
	make grappino;
	make mgp2pdb;
	make mhistogram;
	make cleanobj;
	make version=MPI;
