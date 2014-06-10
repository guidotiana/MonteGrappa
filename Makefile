CC=gcc
CCMP= mpicc
LFLAGS= -Wall -lgsl -lgslcblas -lm -L/usr/local/lib
CFLAGS= -Wall -g -funroll-all-loops -finline-functions -I/opt/local/include -Iinclude/
SRCDIR=src
OBJDIR=obj
INCL = -I../include
INCMPIFLAG=-D ACTIVE_MPI


ifneq ($(version),MPI)


$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< $(INCL) -o $@
	


montegrappa:  $(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/mc.o $(OBJDIR)/local_move.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/potential.o $(OBJDIR)/optimizepot.o $(OBJDIR)/montegrappa.o 
	$(CC) $(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/mc.o  $(OBJDIR)/local_move.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/potential.o $(OBJDIR)/optimizepot.o $(OBJDIR)/montegrappa.o -o montegrappa $(LFLAGS)
	

else

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CCMP) $(INCMPIFLAG) $(CFLAGS) -c $< $(INCL) -o $@




montegrappa_MPI:    $(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/mc.o $(OBJDIR)/local_move.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/potential.o $(OBJDIR)/optimizepot.o $(OBJDIR)/montegrappa.o $(OBJDIR)/MPIfunc.o
	$(CCMP) $(INCMPIFLAG) $(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/mc.o $(OBJDIR)/local_move.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/potential.o $(OBJDIR)/optimizepot.o $(OBJDIR)/montegrappa.o $(OBJDIR)/MPIfunc.o -o montegrappa $(LFLAGS) 

endif

grappino:    $(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/potential.o $(OBJDIR)/energy.o $(OBJDIR)/pdb.o $(OBJDIR)/rotamers.o $(OBJDIR)/grappino.o $(OBJDIR)/optimizepot.o 
	$(CC) $(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/potential.o $(OBJDIR)/energy.o $(OBJDIR)/pdb.o $(OBJDIR)/rotamers.o $(OBJDIR)/optimizepot.o $(OBJDIR)/grappino.o -o grappino $(LFLAGS)


clean:
	rm -f $(OBJDIR)/*.o montegrappa grappino



