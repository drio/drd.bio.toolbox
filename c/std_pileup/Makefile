CC=gcc
CFLAGS=-g -Wall -O2
DFLAGS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE
SAMTOOLS=$(SAMTOOLS_DIR)
INCLUDES=-I$(SAMTOOLS)
TOOL=std_pileup
LIBPATH=-L$(SAMTOOLS)
BIN="../bin"

all: std_pileup install

std_pileup: std_pileup.c
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDES) $(LIBPATH) -lbam -lm -lz $< -lbam -o $@

install:
	cp $(TOOL) $(BIN)

clean:
	rm -f *.o $(TOOL)	
