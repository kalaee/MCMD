# C definitions
CC = gcc
CFLAGS = -Wall -lm -std=gnu99 -O2
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

# abbreviations
CORE = /nfs/users2/kalaee/Documents/core
CORE_= ${CORE}/core_
COULOMB=${CORE}/coulomb_matrix

# object files
${CORE_}rstat.o: ${CORE_}rstat.c ${CORE_}rstat.h
	$(CC) $(CFLAGS) -c $< -o $@

# routines
DEPEND = ${CORE_}rstat.o
parton: parton.c $(DEPEND)
	$(CC) $(DEPEND) $< $(CFLAGS) $(LDLIBS) -o $@
clean:
	rm -f *.o
