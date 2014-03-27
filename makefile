CC=mpicc
CFLAGS= -O2 -Wall -lm


all: ga

ga: fitness.o ga.o init.o mt_mpi.o report.o
    $(CC) $(CFLAGS) fitness.o ga.o init.o mt_mpi.o report.o -o ga
ga.o: ga.c config.h fitness.h ga.h init.h mt_mpi.h report.h types.h
    $(CC) $(CFLAGS) -c ga.c

fitness.o: fitness.c config.h fitness.h types.h
    $(CC) $(CFLAGS) -c fitness.c

init.o: init.c config.h init.h types.h mt_mpi.h report.h
    $(CC) $(CFLAGS) -c init.c
    
mt_mpi.o: mt_mpi.c mt_mpi.h
    $(CC) $(CFLAGS) -c mt_mpi.c
    
report.o: report.c config.h report.h types.h
    $(CC) $(CFLAGS) -c report.c
    

.PSEUDO: clean distclean

clean:
    rm *.o
distclean:
    rm *.o ga
