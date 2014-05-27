CC= mpicc
LFLAGS= -Wall -lgsl -lgslcblas -lm -L/usr/local/lib
CFLAGS= -Wall -g -funroll-all-loops -finline-functions -I/opt/local/include 

.o:     $(CC)  $(CFLAGS) $<

montegrappa:	geometry.o io.o mc.o bias_mp.o memory.o misc.o potential.o optimizepot.o montegrappa.o stempering/stempering.o stempering/memory.o stempering/adjust_st.o stempering/do_mhistogram.o MPIfunc.o
	$(CC) geometry.o io.o mc.o bias_mp.o memory.o misc.o potential.o optimizepot.o montegrappa.o stempering/stempering.o stempering/memory.o stempering/adjust_st.o stempering/do_mhistogram.o MPIfunc.o\
		-o montegrappa $(LFLAGS)

clean:
	rm -f *.o montegrappa grappino stempering/*.o

grappino:	geometry.o io.o  memory.o misc.o potential.o energy.o pdb.o rotamers.o grappino.o optimizepot.o MPIfunc.o mc.o stempering/stempering.o stempering/memory.o stempering/adjust_st.o stempering/do_mhistogram.o
	$(CC) geometry.o io.o  memory.o misc.o potential.o energy.o pdb.o rotamers.o grappino.o optimizepot.o MPIfunc.o mc.o stempering/stempering.o stempering/memory.o stempering/adjust_st.o stempering/do_mhistogram.o -o grappino $(LFLAGS)


