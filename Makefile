CC=gcc
CFLAGS= -Wall -std=gnu99
TUNE= -O 

all: traditional Strassen

traditional:
		$(CC) $(TUNE) $(CFLAGS) -o matrix_mul matrix_mul.c -lpthread

Strassen:
		$(CC) $(TUNE) $(CFLAGS) -o matrix_Strassen matrix_Strassen.c -lpthread