CC= cc
LFLAGS=  -Wall -lgsl -lgslcblas -lm -L/opt/local/lib
CFLAGS=  -Wall -O3 -funroll-all-loops -finline-functions -I/opt/local/include

.o:     $(CC)  $(CFLAGS) $<

clean:
	rm *.o

all: stempering.o memory.o adjust_st.o do_mhistogram.o


